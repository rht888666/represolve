# Minimal runtime for the 30-second demo and package import checks.
FROM python:3.11-slim

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY demo ./demo

RUN pip install --no-cache-dir .

CMD ["python", "demo/run_demo.py"]
