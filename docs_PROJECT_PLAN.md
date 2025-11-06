# BioCLI Companion - Detailed Project Plan

## 🎯 Project Overview

**Goal**: Create a web-based tool that helps students and researchers understand bioinformatics command-line workflows through interactive explanations and guided command building.

**Target Users**: 
- Bioinformatics students (primary)
- Researchers new to computational biology
- Anyone learning bioinformatics pipelines

**Success Metrics**:
- You use it regularly in your own work
- 5+ classmates find it helpful
- Positive feedback on r/bioinformatics
- Portfolio piece for grad school/jobs

---

## 📅 Timeline & Phases

### Phase 1: MVP - Command Explainer (Weeks 1-4)

**Goal**: Build the simplest useful version

#### Week 1: Setup & Prototype
**Deliverable**: Basic working prototype

**Tasks**:
1. ✅ Create GitHub repository
2. [ ] Set up basic project structure
3. [ ] Create simple HTML/JS frontend (single page)
4. [ ] Set up OpenAI API integration
5. [ ] Build command parser (basic)
6. [ ] Test with 3-5 common commands

**Validation**: Can paste `bwa mem -t 8 ref.fa reads.fq` and get a decent explanation

#### Week 2: Backend Foundation & API Issues
**Deliverable**: Python FastAPI backend with working OpenAI integration

**Tasks**:
1. [x] Set up FastAPI project
2. [x] Create API endpoint for command explanation
3. [x] Implement command parsing logic
4. [x] Add bioinformatics tool detection (BWA, samtools, FastQC, etc.)
5. [x] Integrate OpenAI API properly
6. [x] Add error handling
7. [ ] **🔧 FIX OPENAI API QUOTA ISSUE** - Check billing/usage and restore real AI responses
8. [ ] Test real OpenAI responses with multiple commands
9. [ ] Optimize prompt for better bioinformatics explanations

**Validation**: API returns structured explanations reliably

#### Week 3: Frontend Improvement
**Deliverable**: Clean, usable interface

**Tasks**:
1. [ ] Set up React app (or keep it simple with vanilla JS)
2. [ ] Create command input component
3. [ ] Build explanation display component
4. [ ] Add loading states and error handling
5. [ ] Basic styling (Tailwind or simple CSS)
6. [ ] Mobile-responsive design

**Validation**: Interface is intuitive and looks professional

#### Week 4: Polish & Deploy
**Deliverable**: Live, usable tool

**Tasks**:
1. [ ] Add support for 10 common tools
2. [ ] Improve explanation quality (iterate based on testing)
3. [ ] Add example commands
4. [ ] Write basic documentation
5. [ ] Deploy frontend (Vercel)
6. [ ] Deploy backend (Railway or Render)
7. [ ] Test with 5 users

**Validation**: Tool is live and receiving positive feedback

---

### Phase 2: Command Builder (Weeks 5-8)

**Goal**: Add interactive command generation

#### Week 5: Design & Planning
**Tasks**:
1. [ ] Design command builder UX
2. [ ] Map out parameter forms for each tool
3. [ ] Create tool configuration files
4. [ ] Plan database schema for storing commands

#### Week 6-7: Implementation
**Tasks**:
1. [ ] Build command builder UI
2. [ ] Create parameter input forms (per tool)
3. [ ] Implement command generation logic
4. [ ] Add validation for parameters
5. [ ] Create "Export to script" feature
6. [ ] Add command history

#### Week 8: Integration & Testing
**Tasks**:
1. [ ] Integrate builder with explainer
2. [ ] Add bi-directional flow (explain → edit → rebuild)
3. [ ] User testing with classmates
4. [ ] Refinement based on feedback

---

### Phase 3: Enhanced Learning (Weeks 9-12)

**Goal**: Add educational features

#### Features to Add:
1. [ ] Common error detection
   - Parse error messages
   - Suggest solutions
   - Link to relevant documentation

2. [ ] Best practices recommendations
   - "Consider adding -R for read groups"
   - "This will use X GB of RAM"
   - "Tip: Index your reference first"

3. [ ] Interactive tutorials
   - Step-by-step workflows (RNA-seq, variant calling, etc.)
   - Progress tracking
   - Quiz/validation questions

4. [ ] Resource integration
   - Link to official documentation
   - Video tutorials
   - Community forum discussions

---

## 🛠️ Technical Architecture

### System Components

```
┌─────────────────────────────────────────────┐
│              Frontend (React)                │
│  ┌────────────┐         ┌─────────────────┐ │
│  │  Explainer │         │ Command Builder │ │
│  └────────────┘         └─────────────────┘ │
└─────────────────┬───────────────────────────┘
                  │ REST API
┌─────────────────▼───────────────────────────┐
│          Backend (FastAPI)                   │
│  ┌──────────────┐  ┌────────────────────┐   │
│  │Command Parser│  │ AI Explanation Gen │   │
│  └──────────────┘  └────────────────────┘   │
│  ┌──────────────┐  ┌────────────────────┐   │
│  │Tool Validator│  │  Error Detector    │   │
│  └──────────────┘  └────────────────────┘   │
└─────────────────┬───────────────────────────┘
                  │
        ┌─────────┴─────────┐
        │                   │
┌───────▼────────┐  ┌──────▼──────┐
│ OpenAI GPT-4   │  │   Database  │
│      API       │  │   (SQLite)  │
└────────────────┘  └─────────────┘
```

### API Endpoints (Phase 1)

```
POST /api/explain
- Input: { "command": "string" }
- Output: { 
    "tool": "bwa mem",
    "explanation": {...},
    "warnings": [...],
    "resources": [...]
  }

GET /api/tools
- Output: List of supported tools

GET /api/examples
- Output: Example commands for learning
```

### Database Schema (Phase 2)

```sql
-- User sessions (optional, for tracking progress)
CREATE TABLE sessions (
    id INTEGER PRIMARY KEY,
    session_id TEXT UNIQUE,
    created_at TIMESTAMP
);

-- Command history
CREATE TABLE command_history (
    id INTEGER PRIMARY KEY,
    session_id TEXT,
    command TEXT,
    timestamp TIMESTAMP,
    FOREIGN KEY (session_id) REFERENCES sessions(session_id)
);

-- Tool configurations
CREATE TABLE tools (
    id INTEGER PRIMARY KEY,
    name TEXT UNIQUE,
    description TEXT,
    parameters JSON,
    examples JSON
);
```

---

## 📚 Supported Tools Roadmap

### Phase 1 (MVP): Core 10 Tools
1. **bwa** - Read alignment
2. **samtools** - SAM/BAM manipulation
3. **FastQC** - Quality control
4. **STAR** - RNA-seq alignment
5. **bcftools** - Variant calling
6. **bedtools** - Genomic interval operations
7. **featureCounts** - Read counting
8. **trimmomatic** - Read trimming
9. **picard** - BAM processing
10. **GATK** - Variant discovery

### Phase 2: Extended Tools (15 more)
- bowtie2, hisat2, kallisto, salmon, cutadapt, etc.

### Phase 3: Advanced Pipelines
- Full workflow explanations (end-to-end RNA-seq, variant calling, etc.)

---

## 🎨 Design Principles

1. **Clarity over completeness** - Better to explain 10 tools well than 100 poorly
2. **Learning-focused** - Teach, don't just automate
3. **Show, don't hide** - Always show the actual command, don't abstract it away
4. **Progressive disclosure** - Simple by default, detailed on demand
5. **Mobile-friendly** - Should work on phones/tablets for on-the-go learning

---

## 🚀 Deployment Strategy

### Phase 1 (MVP)
- **Frontend**: Vercel (free tier)
  - Automatic deployments from GitHub
  - Global CDN
  - Custom domain support

- **Backend**: Railway or Render (free tier)
  - Easy Python deployment
  - Environment variable management
  - Auto-scaling

- **Database**: SQLite file (simple start)
  - Stored in backend container
  - Easy to backup
  - Sufficient for early users

### Phase 2 (Growth)
- Upgrade to PostgreSQL (Railway/Render managed)
- Consider caching layer (Redis)
- Monitor API usage and costs

---

## 📊 Success Metrics & Validation

### Phase 1 Validation (Week 4)
- [ ] 5+ people use it and find it helpful
- [ ] Positive feedback on usability
- [ ] Explanations are accurate and clear
- [ ] Tool handles 10 common commands reliably

### Phase 2 Validation (Week 8)
- [ ] 20+ active users
- [ ] Users generate commands successfully
- [ ] 80%+ accuracy on generated commands
- [ ] Feedback indicates it speeds up learning

### Phase 3 Validation (Week 12)
- [ ] 50+ users
- [ ] Shared on r/bioinformatics with positive reception
- [ ] Used in at least one workshop/class
- [ ] Clear evidence it helps people learn CLI

---

## 🔄 Iteration Process

**Weekly Cycle**:
1. **Monday**: Set goals for the week
2. **Wed-Fri**: Build features
3. **Weekend**: Test with users, gather feedback
4. **Sunday**: Review what worked, plan next week

**User Testing Checklist**:
- [ ] Watch someone use it (don't help)
- [ ] Ask: "What's confusing?"
- [ ] Ask: "What would make this more useful?"
- [ ] Note where they struggle
- [ ] Prioritize fixes based on feedback

---

## 🎓 Learning Opportunities

This project teaches you:
- Full-stack web development
- API design and implementation
- AI integration (GPT-4)
- User experience design
- DevOps and deployment
- Community building
- Product iteration

**Bonus**: Portfolio piece for grad school applications and job interviews!

---

## 🔮 Future Possibilities

**If this gains traction**, consider:
- VS Code extension version
- Slack/Discord bot for bioinformatics communities
- Integration with Jupyter notebooks
- Workflow validation service
- Commercial features (team collaboration, custom tools)

**But focus on MVP first!** Don't get distracted by shiny features.

---

## 📝 Notes & Decisions

### Technology Choices Explained

**Why React?**
- Industry standard, good for portfolio
- Component reusability
- Large community for help
- (Alternative: Keep it simple with vanilla JS for MVP)

**Why FastAPI?**
- Fast to develop (Python)
- Automatic API documentation
- Async support for AI API calls
- Type hints improve code quality

**Why OpenAI?**
- Best explanations quality
- Easy to integrate
- Cost is manageable for this use case ($5-20/month in testing)
- Can switch to open models later if needed

---

## ✅ Current Status

**Phase**: Phase 1 - Week 1
**Next Milestone**: Working prototype (end of Week 1)
**Blockers**: None currently

---

**Remember**: The goal is to help people understand CLI, not to build the perfect tool. Ship early, iterate based on feedback, and focus on making it genuinely useful for learners.