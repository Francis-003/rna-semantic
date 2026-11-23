"""可选 LLM 润色：仅用标准库调用 OpenAI 兼容 Chat Completions API。

安全约束：
- 输入仅为结构化事实 JSON，不含原始数据；
- system prompt 强制"只润色语言、保留全部数值与 P/Z 锚点、不得新增事实或临床结论"；
- 无 API key 或调用失败时自动降级为规则版文本，管线不中断。
- SSRF 防护：仅允许 HTTPS + 域名白名单；解析出的所有 IP 必须为公网；
  连接时固定到已校验 IP（防 DNS rebinding），SNI/证书校验仍用原域名；
  禁用代理与重定向。
"""
from __future__ import annotations

import http.client
import ipaddress
import json
import socket
import ssl
import urllib.error
import urllib.parse
import urllib.request

from . import config as cfg

SYSTEM_PROMPT = """你是一名谨慎的生物医学文本编辑。你的唯一任务是润色一段已有的结构化临床摘要，使其语言更流畅、专业、符合中文临床表达习惯。

严格规则（违反即视为失败）：
1. 只能使用输入中已提供的事实、数值、百分位(Pxx)和 Z-score，一个都不能删改，一个都不能新增。
2. 不得添加任何输入中没有的基因、数值、结论、因果关系或临床建议。
3. 不得把 mRNA 表达推断为 DNA 突变；涉及扩增只能沿用原文的"提示可能存在扩增/过表达"审慎措辞。
4. 必须完整保留五个分节标题与全部量化锚点（如 P96、Z=+2.31、MT%=8.2%）。
5. 保留"临床解读边界"段落的免责声明含义，不得削弱。
6. 不要改变 Markdown 结构（标题层级、列表符号、粗体标记）。
只输出润色后的 Markdown 文本，不要任何额外解释。"""

# SSRF 防护：仅允许下列 HTTPS 主机（OpenAI 官方及常见兼容网关域名）
ALLOWED_HOSTS = {
    "api.openai.com",
    "api.anthropic.com",
    "openai.googleapis.com",
}

_PRIVATE_NETWORKS = [
    ipaddress.ip_network(n) for n in (
        "0.0.0.0/8", "10.0.0.0/8", "100.64.0.0/10", "127.0.0.0/8",
        "169.254.0.0/16", "172.16.0.0/12", "192.0.0.0/24", "192.168.0.0/16",
        "198.18.0.0/15", "224.0.0.0/4", "240.0.0.0/4", "255.255.255.255/32",
        "::1/128", "fc00::/7", "fe80::/10", "::ffff:0:0/96",
        "2001:db8::/32",  # 文档保留地址
    )
]


def _is_public_ip(ip: str) -> bool:
    try:
        addr = ipaddress.ip_address(ip)
    except ValueError:
        return False
    if addr.is_unspecified or addr.is_loopback or addr.is_link_local or addr.is_reserved:
        return False
    if addr.is_private:
        return False
    return not any(addr in net for net in _PRIVATE_NETWORKS)


def _validate_url(base_url: str):
    """返回 (api_url, host, pinned_ip)；不安全则返回 None。"""
    try:
        parsed = urllib.parse.urlparse(base_url)
    except ValueError:
        return None
    if parsed.scheme != "https":
        return None
    host = (parsed.hostname or "").lower()
    if host not in ALLOWED_HOSTS:
        return None
    if parsed.username or parsed.password or parsed.query or parsed.fragment:
        return None
    try:
        infos = socket.getaddrinfo(host, 443, type=socket.SOCK_STREAM)
    except socket.gaierror:
        return None
    ips = {info[4][0] for info in infos}
    if not ips or not all(_is_public_ip(ip) for ip in ips):
        return None
    # 固定到其中一个已校验的 IPv4（优先），避免连接阶段被 DNS rebinding
    ipv4 = sorted(ip for ip in ips if ":" not in ip)
    pinned_ip = ipv4[0] if ipv4 else sorted(ips)[0]
    path = parsed.path.rstrip("/")
    if path.endswith("/chat/completions"):
        api_url = f"https://{host}{path}"
    else:
        api_url = f"https://{host}{path}/chat/completions"
    return api_url, host, pinned_ip


def _build_opener(host: str, pinned_ip: str):
    """构造仅连接 pinned_ip、禁用代理与重定向的 opener。"""
    class _NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):
            raise urllib.error.HTTPError(req.full_url, code, "拒绝重定向", headers, fp)

    class _PinnedConnection(http.client.HTTPSConnection):
        def connect(self):
            sock = socket.create_connection((pinned_ip, self.port), timeout=self.timeout)
            # TLS SNI 与证书校验使用原域名，而非 IP，保证证书正确验证
            self.sock = self._context.wrap_socket(sock, server_hostname=host)

    class _PinnedHTTPSHandler(urllib.request.HTTPSHandler):
        def https_open(self, req):
            return self.do_open(
                _PinnedConnection, req,
                context=ssl.create_default_context(),
            )

    return urllib.request.build_opener(
        urllib.request.ProxyHandler({}),  # 禁用环境代理，防止绕过 IP 校验
        _NoRedirect(),
        _PinnedHTTPSHandler(),
    )


def polish_markdown(markdown_text: str, facts: dict) -> tuple[str, str]:
    """返回 (文本, 状态)。状态: 'polished' / 'skipped' / 'failed'。"""
    llm = cfg.llm_config()
    if not llm["enabled"]:
        return markdown_text, "skipped"

    validated = _validate_url(llm["base_url"])
    if validated is None:
        print(
            "[llm_polish] LLM_BASE_URL 未通过安全校验（仅允许 HTTPS、公网白名单域名），"
            "降级为规则版文本。"
        )
        return markdown_text, "failed"
    api_url, host, pinned_ip = validated

    payload = {
        "model": llm["model"],
        "temperature": 0.2,
        "messages": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": (
                "【结构化事实（JSON）】\n"
                + json.dumps(facts, ensure_ascii=False)
                + "\n\n【待润色的 Markdown 文本】\n" + markdown_text
            )},
        ],
    }
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {llm['api_key']}",
    }
    req = urllib.request.Request(
        api_url, data=json.dumps(payload).encode("utf-8"),
        headers=headers, method="POST",
    )
    try:
        opener = _build_opener(host, pinned_ip)
        with opener.open(req, timeout=llm["timeout"]) as resp:
            body = json.loads(resp.read().decode("utf-8"))
        text = body["choices"][0]["message"]["content"].strip()
        if not text:
            return markdown_text, "failed"
        return text, "polished"
    except (urllib.error.URLError, OSError, KeyError, IndexError, ValueError) as exc:
        print(f"[llm_polish] 润色失败，降级为规则版文本：{exc}")
        return markdown_text, "failed"
