# Learning Roadmap 📚

This document outlines everything you need to learn to build BioCLI Companion, organized by phase and difficulty.

---

## 🎯 Overview: What You Need to Know

### Already Know (As a Bioinformatics Student)
✅ Python programming  
✅ Bioinformatics tools and workflows  
✅ Command-line basics  
✅ Git/GitHub basics  

### Need to Learn
- Frontend development (React/JavaScript)
- Backend API development (FastAPI)
- AI API integration (OpenAI)
- Deployment and DevOps basics

**Good news**: You can learn as you build!

---

## 📅 Phase 1: MVP (Weeks 1-4)

### Week 1: Quick Start Skills

#### 1. Basic HTML/CSS/JavaScript (4-6 hours)
**Why**: Build a simple frontend quickly

**What to learn**:
- HTML forms and input elements
- CSS for basic styling
- JavaScript DOM manipulation
- Fetch API for HTTP requests

**Resources**:
- [MDN Web Docs - Getting Started](https://developer.mozilla.org/en-US/docs/Learn)
- [JavaScript.info](https://javascript.info/) - First 3 chapters
- [CSS Tricks - Flexbox Guide](https://css-tricks.com/snippets/css/a-guide-to-flexbox/)

**Practice Project**: 
Build a simple calculator web app (1-2 hours)

#### 2. OpenAI API Basics (2-3 hours)
**Why**: Core functionality for explanations

**What to learn**:
- API keys and authentication
- Making API calls to GPT-4
- Prompt engineering basics
- Handling responses

**Resources**:
- [OpenAI API Quickstart](https://platform.openai.com/docs/quickstart)
- [OpenAI Cookbook](https://github.com/openai/openai-cookbook)

**Practice Project**:
Write a Python script that takes a command and gets an explanation from GPT-4 (1 hour)

```python
import openai

openai.api_key = "your-key-here"

def explain_command(command):
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": "You are a bioinformatics CLI expert."},
            {"role": "user", "content": f"Explain this command: {command}"}
        ]
    )
    return response.choices[0].message.content

# Test it
print(explain_command("bwa mem -t 8 ref.fa reads.fq"))
```

#### 3. Git/GitHub Workflow (1-2 hours refresher)
**Why**: Version control and collaboration

**What to learn**:
- Creating branches
- Committing changes
- Pushing to GitHub
- Basic PR workflow

**Resources**:
- [GitHub Skills](https://skills.github.com/)
- [Git Cheat Sheet](https://education.github.com/git-cheat-sheet-education.pdf)

---

### Week 2: Backend Development

#### 1. FastAPI Fundamentals (6-8 hours)
**Why**: Build the backend API

**What to learn**:
- FastAPI project structure
- Creating API endpoints (GET, POST)
- Request/response models (Pydantic)
- CORS configuration
- Error handling
- Environment variables (.env files)

**Resources**:
- [FastAPI Tutorial](https://fastapi.tiangolo.com/tutorial/) - Complete the full tutorial
- [Real Python - FastAPI Guide](https://realpython.com/fastapi-python-web-apis/)

**Practice Project**: Build a simple note-taking API (3-4 hours)
```python
# main.py
from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI()

class Note(BaseModel):
    title: str
    content: str

notes = []

@app.post("/notes")
def create_note(note: Note):
    notes.append(note)
    return {"status": "created", "note": note}

@app.get("/notes")
def get_notes():
    return {"notes": notes}
```

#### 2. Python Virtual Environments (1 hour)
**Why**: Dependency management

**What to learn**:
- Creating virtual environments
- Installing packages with pip
- requirements.txt files

**Resources**:
- [Python venv documentation](https://docs.python.org/3/library/venv.html)

**Practice**:
```bash
python -m venv venv
source venv/bin/activate
pip install fastapi uvicorn openai
pip freeze > requirements.txt
```

#### 3. API Testing (2-3 hours)
**Why**: Ensure your API works correctly

**What to learn**:
- Using Postman or Thunder Client
- Writing basic pytest tests
- FastAPI's automatic /docs endpoint

**Resources**:
- [FastAPI Testing](https://fastapi.tiangolo.com/tutorial/testing/)
- [Postman Learning Center](https://learning.postman.com/)

---

### Week 3: Frontend Development

#### 1. React Basics (8-10 hours)
**Why**: Build a modern, interactive UI

**What to learn**:
- Components and props
- State management (useState)
- Effects (useEffect)
- Event handling
- Conditional rendering
- Making API calls from React

**Resources**:
- [React Official Tutorial](https://react.dev/learn) - Tic-Tac-Toe tutorial
- [React Beta Docs](https://react.dev/) - Read "Learn React" section
- [Scrimba React Course](https://scrimba.com/learn/learnreact) - Free, interactive

**Practice Project**: Build a weather app that fetches data from an API (4-5 hours)

**Alternative** (if React feels overwhelming):
Start with vanilla JavaScript and upgrade later. React is nice-to-have, not required for MVP.

#### 2. Tailwind CSS (2-3 hours)
**Why**: Quick, professional styling

**What to learn**:
- Utility-first CSS concepts
- Common Tailwind classes
- Responsive design
- Dark mode

**Resources**:
- [Tailwind CSS Docs](https://tailwindcss.com/docs)
- [Tailwind UI Components](https://tailwindui.com/components) (free examples)

**Practice**: 
Recreate a simple design from Dribbble using Tailwind (2 hours)

#### 3. Modern JavaScript (ES6+) (3-4 hours)
**Why**: Write cleaner React code

**What to learn**:
- Arrow functions
- Destructuring
- Template literals
- Async/await
- Array methods (map, filter, reduce)
- Modules (import/export)

**Resources**:
- [JavaScript.info - Modern JavaScript](https://javascript.info/)
- [ES6 Features](http://es6-features.org/)

---

### Week 4: Deployment

#### 1. Deploying FastAPI (3-4 hours)
**Why**: Get your backend online

**What to learn**:
- Environment variables in production
- Railway or Render deployment
- Monitoring logs
- HTTPS/SSL basics

**Resources**:
- [Deploying FastAPI to Railway](https://railway.app/template/fastapi)
- [Render Python Deployment](https://render.com/docs/deploy-fastapi)

**Steps**:
1. Create account on Railway/Render
2. Connect GitHub repo
3. Set environment variables (OpenAI key)
4. Deploy and test

#### 2. Deploying React to Vercel (2-3 hours)
**Why**: Get your frontend online

**What to learn**:
- Build process (npm run build)
- Environment variables in Vercel
- Custom domains (optional)
- Connecting to backend API

**Resources**:
- [Vercel Deployment Docs](https://vercel.com/docs)
- [Deploy React to Vercel](https://vercel.com/guides/deploying-react-with-vercel)

**Steps**:
1. Create Vercel account
2. Import GitHub repo
3. Configure build settings
4. Deploy

#### 3. CORS Understanding (1-2 hours)
**Why**: Allow frontend to talk to backend

**What to learn**:
- What CORS is and why it exists
- Configuring CORS in FastAPI
- Allowed origins in production

**Resources**:
- [MDN - CORS](https://developer.mozilla.org/en-US/docs/Web/HTTP/CORS)
- [FastAPI CORS](https://fastapi.tiangolo.com/tutorial/cors/)

---

## 📅 Phase 2: Command Builder (Weeks 5-8)

### Advanced Frontend Skills

#### 1. Forms and Validation (4-5 hours)
**What to learn**:
- Controlled components
- Form libraries (React Hook Form)
- Input validation
- Error messages

**Resources**:
- [React Hook Form](https://react-hook-form.com/)
- [Form Validation Tutorial](https://www.freecodecamp.org/news/how-to-validate-forms-in-react/)

#### 2. State Management (3-4 hours)
**What to learn**:
- Context API (for global state)
- When to use local vs. global state
- State organization patterns

**Resources**:
- [React Context API](https://react.dev/learn/passing-data-deeply-with-context)

#### 3. Dynamic UI Components (3-4 hours)
**What to learn**:
- Rendering forms dynamically based on tool
- Component composition
- Reusable input components

---

### Database Basics

#### 1. SQLite with Python (4-5 hours)
**What to learn**:
- SQL basics (SELECT, INSERT, UPDATE, DELETE)
- Using SQLite with Python (sqlite3 module)
- SQLAlchemy ORM basics

**Resources**:
- [SQLite Tutorial](https://www.sqlitetutorial.net/)
- [SQLAlchemy Tutorial](https://docs.sqlalchemy.org/en/20/tutorial/)

**Practice**:
Build a simple TODO app with database persistence (3 hours)

---

## 📅 Phase 3: Enhanced Features (Weeks 9-12)

### Advanced AI Integration

#### 1. Prompt Engineering (5-6 hours)
**What to learn**:
- Writing effective prompts
- Few-shot learning
- Chain-of-thought prompting
- Handling edge cases

**Resources**:
- [OpenAI Prompt Engineering Guide](https://platform.openai.com/docs/guides/prompt-engineering)
- [Prompt Engineering Guide](https://www.promptingguide.ai/)

#### 2. Error Handling & Parsing (4-5 hours)
**What to learn**:
- Regex for parsing error messages
- Pattern matching
- Structured error responses

---

## 🎓 Optional Advanced Topics

### If You Want to Go Deeper

#### 1. TypeScript (10-15 hours)
**Why**: Better code quality and IDE support

**When**: After you're comfortable with JavaScript

**Resources**:
- [TypeScript Handbook](https://www.typescriptlang.org/docs/handbook/intro.html)

#### 2. Testing (8-10 hours)
**Why**: Ensure code reliability

**What to learn**:
- pytest for backend
- React Testing Library for frontend
- Test-driven development basics

**Resources**:
- [pytest Documentation](https://docs.pytest.org/)
- [React Testing Library](https://testing-library.com/docs/react-testing-library/intro/)

#### 3. CI/CD (4-6 hours)
**Why**: Automated testing and deployment

**What to learn**:
- GitHub Actions basics
- Automated testing on commits
- Automated deployments

**Resources**:
- [GitHub Actions Documentation](https://docs.github.com/en/actions)

---

## 📚 General Learning Resources

### Best Practices
- [Clean Code Principles](https://github.com/ryanmcdermott/clean-code-javascript)
- [Python Best Practices](https://realpython.com/tutorials/best-practices/)
- [React Patterns](https://reactpatterns.com/)

### Communities for Help
- [r/bioinformatics](https://www.reddit.com/r/bioinformatics/)
- [r/learnprogramming](https://www.reddit.com/r/learnprogramming/)
- [Stack Overflow](https://stackoverflow.com/)
- [Discord - Reactiflux](https://www.reactiflux.com/) (React help)
- [Discord - Python](https://pythondiscord.com/)

### Documentation to Bookmark
- [MDN Web Docs](https://developer.mozilla.org/)
- [React Docs](https://react.dev/)
- [FastAPI Docs](https://fastapi.tiangolo.com/)
- [Python Docs](https://docs.python.org/3/)

---

## 📊 Learning Time Estimates

### Total Time by Phase

**Phase 1 (MVP)**:
- Week 1: ~15 hours learning + 10 hours building
- Week 2: ~12 hours learning + 15 hours building
- Week 3: ~15 hours learning + 15 hours building
- Week 4: ~8 hours learning + 20 hours building

**Total Phase 1**: ~50 hours learning + ~60 hours building = **110 hours** (realistic for one semester)

**Phase 2**: ~30 hours learning + ~50 hours building = **80 hours**

**Phase 3**: ~20 hours learning + ~40 hours building = **60 hours**

---

## 🎯 Learning Strategy Tips

### 1. **Learn by Building**
Don't just read tutorials - build the practice projects. You'll retain more.

### 2. **Focus on Breadth First**
Get a basic understanding of everything before going deep. You can specialize later.

### 3. **Use AI to Help**
- ChatGPT/Claude for explaining concepts
- GitHub Copilot for code suggestions
- Don't copy-paste blindly - understand what the code does

### 4. **Build in Public**
- Share your progress on Twitter/LinkedIn
- Ask questions in communities
- Teaching others reinforces your learning

### 5. **Don't Get Stuck**
- If stuck for >30 minutes, ask for help
- Move on and come back later
- Perfect is the enemy of done (for MVP)

---

## ✅ Weekly Checkpoints

Use this to track your progress:

### Week 1
- [ ] Completed HTML/CSS/JS basics
- [ ] Successfully called OpenAI API
- [ ] Created GitHub repo
- [ ] Built simple prototype

### Week 2
- [ ] Completed FastAPI tutorial
- [ ] Built practice API
- [ ] Integrated OpenAI with FastAPI
- [ ] Tested API with Postman

### Week 3
- [ ] Completed React tutorial
- [ ] Built practice React app
- [ ] Created BioCLI Companion frontend
- [ ] Connected frontend to backend

### Week 4
- [ ] Deployed backend to Railway/Render
- [ ] Deployed frontend to Vercel
- [ ] Fixed CORS issues
- [ ] Tested live site

---

**Remember**: You don't need to be an expert in everything. You just need to know enough to build something useful. Learn as you go, and don't be afraid to ask for help!

**The best way to learn is to build something you'll actually use.**