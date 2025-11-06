# BioCLI Companion Backend

## Quick Start

1. Create virtual environment:
```bash
python -m venv venv
```

2. Activate virtual environment:
```bash
# Windows
venv\Scripts\activate

# macOS/Linux
source venv/bin/activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Set up environment variables:
```bash
copy .env.example .env
# Edit .env file with your OpenAI API key
```

5. Run the server:
```bash
python main.py
```

The API will be available at: http://localhost:8000

## API Documentation

Once running, visit http://localhost:8000/docs for interactive API documentation.

## Environment Variables

- `OPENAI_API_KEY`: Your OpenAI API key
- `API_HOST`: Server host (default: localhost)
- `API_PORT`: Server port (default: 8000)
- `DEBUG`: Enable debug mode (default: True)