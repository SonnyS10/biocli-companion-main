# BioCLI Companion 🧬

> An interactive tool to help students and researchers understand bioinformatics command-line workflows

## 🎯 Mission

Make the command-line portion of bioinformatics easier to understand for people learning computational biology.

## 🚀 What This Does

BioCLI Companion is a web-based tool that helps you:
- **Explain** complex bioinformatics commands in plain English
- **Build** commands step-by-step with guided inputs
- **Understand** what each parameter does and why it matters
- **Learn** by doing, with real-time feedback

## ✨ Features (Planned)

### Phase 1: MVP (Current Focus)
- [ ] Command explainer - paste a command, get detailed breakdown
- [ ] AI-powered explanations using GPT-4
- [ ] Support for common tools (BWA, samtools, FastQC, STAR, etc.)
- [ ] Clean, simple web interface

### Phase 2: Command Builder
- [ ] Interactive form to build commands step-by-step
- [ ] Parameter validation and suggestions
- [ ] Export to shell script

### Phase 3: Enhanced Learning
- [ ] Common error detection and solutions
- [ ] Best practices recommendations
- [ ] Video tutorials and documentation links

## 🛠️ Tech Stack

- **Frontend**: React + Tailwind CSS
- **Backend**: Python (FastAPI)
- **AI**: OpenAI GPT-4 API
- **Database**: SQLite (initially)
- **Deployment**: Vercel (frontend) + Railway/Render (backend)

## 📚 Documentation

- [Project Plan](./docs/PROJECT_PLAN.md) - Detailed roadmap and phases
- [Learning Roadmap](./docs/LEARNING_ROADMAP.md) - What you need to learn for each phase
- [Contributing](./CONTRIBUTING.md) - How to contribute (coming soon)

## 🏃 Quick Start

### Prerequisites
- Python 3.8+
- Node.js 16+
- OpenAI API key

### Installation

```bash
# Clone the repository
git clone https://github.com/SonnyS10/biocli-companion.git
cd biocli-companion

# Set up backend
cd backend
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt

# Set up frontend
cd ../frontend
npm install

# Configure environment variables
cp .env.example .env
# Edit .env and add your OpenAI API key
```

### Running Locally

```bash
# Terminal 1 - Backend
cd backend
source venv/bin/activate
uvicorn main:app --reload

# Terminal 2 - Frontend
cd frontend
npm run dev
```

Visit `http://localhost:3000` to see the app!

## 📖 Example Usage

### Command Explainer
```
Input: bwa mem -t 8 hg38.fa sample_R1.fq sample_R2.fq | samtools sort -o output.bam

Output:
🧬 This command performs paired-end read alignment and sorting:

1. bwa mem - BWA-MEM alignment algorithm
   ├─ -t 8: Use 8 CPU threads for parallel processing
   ├─ hg38.fa: Reference genome (human genome build 38)
   ├─ sample_R1.fq: Forward reads
   └─ sample_R2.fq: Reverse reads

2. | - Pipe output directly to next command (saves disk space)

3. samtools sort - Sort alignments by genomic position
   └─ -o output.bam: Write sorted output to this file

💡 Why this matters:
- Sorted BAM files are required for downstream analysis
- Piping avoids creating large intermediate SAM files
- Using threads (-t 8) speeds up alignment significantly

⚠️ Common issues:
- Ensure reference genome is indexed (bwa index hg38.fa)
- Check you have enough RAM (rule of thumb: 2GB per thread)
- Verify input files are valid FASTQ format
```

## 🎓 Target Audience

- Bioinformatics students learning command-line tools
- Wet lab researchers transitioning to computational analysis
- Anyone who finds bioinformatics CLI intimidating

## 🤝 Contributing

This is an open-source educational project! Contributions are welcome.

Areas where help is needed:
- Adding support for more bioinformatics tools
- Improving explanations and error messages
- UI/UX improvements
- Documentation and tutorials

## 📝 License

MIT License - see [LICENSE](./LICENSE) for details

## 🙏 Acknowledgments

Built by [@SonnyS10](https://github.com/SonnyS10) as a tool to make bioinformatics more accessible.

Inspired by the struggle of learning command-line tools and the desire to help others overcome that barrier.

## 📬 Contact

Questions? Suggestions? Open an issue or reach out!

---

**Status**: 🚧 Under active development - MVP coming soon!