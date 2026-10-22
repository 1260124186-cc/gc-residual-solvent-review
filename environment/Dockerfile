FROM python:3.12.8-slim-bookworm@sha256:2199a62885a12290dc9c5be3ca0681d367576ab7bf037da120e564723292a2f0

RUN python - <<'PY'
import hashlib
import io
import os
import platform
import tarfile
import urllib.request

version = "22.23.4"
arch_map = {
    "x86_64": ("x64", "1805bcd35e6815559bd2eb97db438d3fdb0eb3a5e6a977b7c1badca76337e9ce"),
    "aarch64": ("arm64", "6287aa46b6c375d816e63ae66d024fbbfd70404825aabdfdf898a970871d2f4b"),
}
arch, expected = arch_map[platform.machine()]
filename = f"node-v{version}-linux-{arch}.tar.xz"
url = f"https://npmmirror.com/mirrors/node/v{version}/{filename}"
with urllib.request.urlopen(url) as response:
    archive_bytes = response.read()
if hashlib.sha256(archive_bytes).hexdigest() != expected:
    raise SystemExit("Node.js archive checksum mismatch")
with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:xz") as archive:
    archive.extractall("/opt", filter="data")
os.rename(f"/opt/node-v{version}-linux-{arch}", "/opt/node")
PY

ENV PATH="/opt/node/bin:${PATH}"

RUN npm install --registry=https://registry.npmmirror.com -g @openai/codex@0.163.0 \
    && rm -rf /root/.npm \
    && ln -sf /opt/node/bin/node /usr/local/bin/node \
    && ln -sf /opt/node/bin/npm /usr/local/bin/npm \
    && ln -sf /opt/node/bin/codex /usr/local/bin/codex

RUN find /usr/local/lib/python3.12 -type f -name '*.pyc' -delete \
    && find /usr/local/lib/python3.12 -type d -name '__pycache__' -empty -delete

RUN useradd --create-home --uid 10001 worker
WORKDIR /workspace
COPY workspace/ /workspace/
RUN chown -R worker:worker /workspace
USER worker

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

CMD ["bash"]
