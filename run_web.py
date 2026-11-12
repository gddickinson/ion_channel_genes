"""Start the read-only results dashboard at http://127.0.0.1:8180 (local only).

Run: python run_web.py   (needs requirements-web.txt)
"""
import uvicorn

if __name__ == "__main__":
    uvicorn.run("web.app:app", host="127.0.0.1", port=8180)
