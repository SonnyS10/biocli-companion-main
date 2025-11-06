# Getting Started 🚀

Quick guide to get BioCLI Companion running on your machine.

## Prerequisites

Before you start, make sure you have:
- **Python 3.8+** installed ([Download](https://www.python.org/downloads/))
- **Node.js 16+** installed ([Download](https://nodejs.org/))
- **Git** installed ([Download](https://git-scm.com/))
- **OpenAI API key** ([Get one here](https://platform.openai.com/api-keys))

### Check Your Versions

```bash
python --version  # Should be 3.8 or higher
node --version    # Should be 16 or higher
npm --version
git --version
```

---

## Step 1: Clone the Repository

```bash
git clone https://github.com/SonnyS10/biocli-companion.git
cd biocli-companion
```

---

## Step 2: Set Up the Backend

### Create Virtual Environment

```bash
cd backend
python -m venv venv
```

### Activate Virtual Environment

**On macOS/Linux:**
```bash
source venv/bin/activate
```

**On Windows:**
```bash
venv\Scripts\activate
```

You should see `(venv)` in your terminal prompt.

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Configure Environment Variables

Create a `.env` file in the `backend` directory:

```bash
# On macOS/Linux
touch .env

# On Windows
type nul > .env
```

Add your OpenAI API key to `.env`:
```
OPENAI_API_KEY=your-api-key-here
```

**⚠️ Never commit this file to Git!** It's already in `.gitignore`.

### Test the Backend

```bash
uvicorn main:app --reload
```

Visit http://localhost:8000/docs - you should see the FastAPI documentation.

Press `Ctrl+C` to stop the server.

---

## Step 3: Set Up the Frontend

Open a **new terminal window** (keep the backend running in the other one).

```bash
cd frontend
npm install
```

### Configure Frontend Environment

Create a `.env` file in the `frontend` directory:

```bash
# On macOS/Linux
touch .env

# On Windows
type nul > .env
```

Add your backend URL:
```
REACT_APP_API_URL=http://localhost:8000
```

### Test the Frontend

```bash
npm start
```

The app should open automatically at http://localhost:3000

---

## Step 4: Verify Everything Works

1. **Backend is running**: http://localhost:8000/docs shows API docs
2. **Frontend is running**: http://localhost:3000 shows the app
3. **They can talk**: Try explaining a command in the app

### Test Command

Try pasting this command in the app:
```
bwa mem -t 8 hg38.fa reads_R1.fq reads_R2.fq
```

You should get a detailed explanation!

---

## Common Issues & Solutions

### Issue: "Module not found" errors

**Solution**: Make sure you activated the virtual environment
```bash
source venv/bin/activate  # macOS/Linux
venv\Scripts\activate     # Windows
```

### Issue: "Port already in use"

**Solution**: Kill the process using that port
```bash
# On macOS/Linux
lsof -ti:8000 | xargs kill

# On Windows
netstat -ano | findstr :8000
taskkill /PID <PID> /F
```

### Issue: CORS errors in browser

**Solution**: Make sure backend is configured for your frontend URL in `main.py`:
```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### Issue: OpenAI API errors

**Solutions**:
1. Check your API key is correct in `.env`
2. Verify you have credits: https://platform.openai.com/usage
3. Check API status: https://status.openai.com/

---

## Development Workflow

### Daily Development

**Terminal 1 - Backend:**
```bash
cd backend
source venv/bin/activate  # or venv\Scripts\activate on Windows
uvicorn main:app --reload
```

**Terminal 2 - Frontend:**
```bash
cd frontend
npm start
```

### Making Changes

1. **Backend changes**: Save the file, FastAPI auto-reloads
2. **Frontend changes**: Save the file, React auto-reloads
3. **New dependencies**:
   - Backend: `pip install package-name` then `pip freeze > requirements.txt`
   - Frontend: `npm install package-name`

### Git Workflow

```bash
# Create a new branch for your feature
git checkout -b feature/awesome-feature

# Make changes, then commit
git add .
git commit -m "Add awesome feature"

# Push to GitHub
git push origin feature/awesome-feature
```

---

## Next Steps

1. **Read the docs**:
   - [Project Plan](./PROJECT_PLAN.md) - What we're building
   - [Learning Roadmap](./LEARNING_ROADMAP.md) - What you need to learn

2. **Start building**:
   - Follow Week 1 tasks in the Project Plan
   - Build the practice projects in the Learning Roadmap

3. **Join the community**:
   - Star the repo ⭐
   - Share your progress
   - Ask questions in Issues

---

## Useful Commands Reference

### Backend
```bash
# Activate environment
source venv/bin/activate

# Run server
uvicorn main:app --reload

# Run tests (when you add them)
pytest

# Install new package
pip install package-name
pip freeze > requirements.txt
```

### Frontend
```bash
# Run development server
npm start

# Build for production
npm run build

# Run tests
npm test

# Install new package
npm install package-name
```

### Git
```bash
# Check status
git status

# Create branch
git checkout -b branch-name

# Stage changes
git add .

# Commit
git commit -m "message"

# Push
git push origin branch-name

# Pull latest
git pull origin main
```

---

## Getting Help

- **Stuck?** Open an issue: https://github.com/SonnyS10/biocli-companion/issues
- **Questions?** Check existing issues first
- **Want to contribute?** See CONTRIBUTING.md (coming soon)

---

Happy coding! 🧬