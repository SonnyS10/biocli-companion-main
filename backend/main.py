"""
BioCLI Companion Backend
FastAPI server for explaining bioinformatics commands
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import os
from dotenv import load_dotenv
from openai import OpenAI

# Load environment variables
load_dotenv()

# Initialize OpenAI client
client = OpenAI(
    api_key=os.getenv("OPENAI_API_KEY")
)

app = FastAPI(
    title="BioCLI Companion API",
    description="API for explaining bioinformatics command-line tools",
    version="0.1.0"
)

# CORS middleware for frontend communication
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allow all origins for local development
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class CommandRequest(BaseModel):
    """Request model for command explanation"""
    command: str

class CommandExplanation(BaseModel):
    """Response model for command explanation"""
    tool: str
    command: str
    explanation: str
    success: bool
    error_message: str = None

async def generate_explanation(command: str, tool: str) -> str:
    """
    Generate a detailed explanation of a bioinformatics command using OpenAI.
    
    Args:
        command: The full command to explain
        tool: The primary tool being used
        
    Returns:
        HTML-formatted explanation of the command
    """
    prompt = f"""
    You are an expert bioinformatics educator helping students learn command-line tools for genomics and computational biology.
    
    COMMAND TO EXPLAIN: {command}
    TOOL: {tool}
    
    Please provide a comprehensive, educational explanation formatted in clean HTML. Structure your response exactly like this:

    <h3>🔧 Tool: {tool.upper()}</h3>
    <p><strong>Command:</strong> <code>{command}</code></p>
    
    <h4>📖 What does {tool} do?</h4>
    [Explain the biological purpose and computational function]
    
    <h4>⚙️ Parameter Breakdown</h4>
    <ul>
    [List each parameter with biological/computational context]
    </ul>
    
    <h4>📁 Input & Output Files</h4>
    [Explain file formats, what data they contain, and biological significance]
    
    <h4>🧬 Biological Context</h4>
    [Explain where this fits in genomics workflows - RNA-seq, variant calling, quality control, etc.]
    
    <h4>⚠️ Important Notes</h4>
    <ul>
    [Common pitfalls, computational requirements, best practices]
    </ul>
    
    <h4>🔗 Next Steps</h4>
    [What you typically do after this command in a workflow]

    REQUIREMENTS:
    - Use proper HTML tags and structure
    - Include biological context, not just technical details
    - Explain WHY each parameter matters for the biology
    - Mention computational resources (memory, time) when relevant
    - Use beginner-friendly language but stay scientifically accurate
    - Include specific file format details (FASTQ, SAM, BAM, VCF, etc.)
    """
    
    try:
        api_key = os.getenv('OPENAI_API_KEY')
        print(f"🔑 Using API key: {api_key[:20] if api_key else 'None'}...")  # Debug log
        print(f"📝 Command to explain: {command}")
        print(f"🔧 Tool detected: {tool}")
        
        if not api_key:
            raise Exception("No OpenAI API key found in environment variables")
        
        print("🚀 Making OpenAI API call...")
        # Real OpenAI API call - quota restored!
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "You are an expert bioinformatics educator who explains command-line tools clearly and accurately."},
                {"role": "user", "content": prompt}
            ],
            max_tokens=1500,
            temperature=0.7
        )
        print("✅ OpenAI API call successful!")  # Debug log
        return response.choices[0].message.content
        
    except Exception as e:
        print(f"❌ OpenAI API Error: {str(e)}")  # Debug log
        print(f"❌ Error type: {type(e).__name__}")  # Debug log
        # Fallback to mock response when API fails
        return f"""
        <h3>🔧 Tool: {tool}</h3>
        <p><strong>Command:</strong> <code>{command}</code></p>
        <div style="background: #fff3cd; padding: 15px; border-radius: 5px; border-left: 4px solid #ffc107; margin: 15px 0;">
            <strong>⚠️ AI Service Temporarily Unavailable</strong><br>
            We couldn't generate a detailed explanation right now. Error: {str(e)}
        </div>
        <p>The tool <code>{tool}</code> is a bioinformatics command-line utility. Please try again in a moment for a detailed explanation.</p>
        """

def generate_detailed_mock_explanation(command: str, tool: str) -> str:
    """Generate a detailed mock explanation for development/testing"""
    
    explanations = {
        "bwa": {
            "description": "BWA (Burrows-Wheeler Aligner) is a software package for mapping DNA sequences against a large reference genome.",
            "use_case": "Essential for aligning short sequencing reads to a reference genome in genomics workflows.",
            "common_params": {
                "-t": "Number of threads to use for parallel processing",
                "mem": "BWA-MEM algorithm, best for reads longer than 70bp"
            }
        },
        "samtools": {
            "description": "SAMtools provides various utilities for manipulating alignments in the SAM/BAM format.",
            "use_case": "Converting, sorting, indexing, and viewing alignment files.",
            "common_params": {
                "view": "View and convert SAM/BAM files",
                "-b": "Output in BAM format (binary)",
                "-S": "Input is in SAM format"
            }
        },
        "fastqc": {
            "description": "FastQC is a quality control tool for high throughput sequence data.",
            "use_case": "Assessing the quality of raw sequencing data before analysis.",
            "common_params": {
                "-o": "Output directory",
                "--threads": "Number of files to process simultaneously"
            }
        }
    }
    
    tool_info = explanations.get(tool, {
        "description": f"{tool} is a bioinformatics command-line utility.",
        "use_case": "Used in computational biology workflows.",
        "common_params": {}
    })
    
    # Parse parameters from command
    params = command.split()[1:] if len(command.split()) > 1 else []
    
    return f"""
    <h3>🔧 Tool: {tool.upper()}</h3>
    <p><strong>Command:</strong> <code>{command}</code></p>
    
    <div style="background: #e8f5e8; padding: 15px; border-radius: 5px; border-left: 4px solid #28a745; margin: 15px 0;">
        <strong>✅ Mock Response</strong> (OpenAI quota exceeded - using educational mock data)
    </div>
    
    <h4>📖 What does {tool} do?</h4>
    <p>{tool_info['description']}</p>
    
    <h4>🔬 Common Use Case</h4>
    <p>{tool_info['use_case']}</p>
    
    <h4>⚙️ Parameters in Your Command</h4>
    <ul>
        {''.join([f"<li><code>{param}</code> - Parameter analysis (detailed explanation would come from AI)</li>" for param in params[:3]])}
    </ul>
    
    <h4>📁 Input/Output Files</h4>
    <p>Based on the command structure, this appears to work with sequence data files. Detailed file format explanations would be provided by the AI service.</p>
    
    <h4>⚠️ Important Notes</h4>
    <ul>
        <li>This is a mock response for development purposes</li>
        <li>Real explanations will be much more detailed once OpenAI quota is available</li>
        <li>Each parameter would be explained in biological context</li>
        <li>Workflow integration and best practices would be included</li>
    </ul>
    
    <div style="background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 15px;">
        <small><strong>💡 Next Steps:</strong> Add billing to your OpenAI account or wait for quota reset to get full AI-powered explanations!</small>
    </div>
    """

@app.get("/")
async def root():
    """Health check endpoint"""
    return {"message": "BioCLI Companion API is running!", "version": "0.1.0"}

@app.post("/api/explain", response_model=CommandExplanation)
async def explain_command(request: CommandRequest):
    """
    Explain a bioinformatics command in detail.
    
    Args:
        request: CommandRequest containing the command to explain
        
    Returns:
        Detailed explanation including tool, parameters, and usage
    """
    try:
        command = request.command.strip()
        
        if not command:
            raise HTTPException(status_code=400, detail="Command cannot be empty")
        
        # Extract the main tool from the command
        tool = command.split()[0] if command.split() else "unknown"
        
        # Generate explanation using OpenAI
        explanation = await generate_explanation(command, tool)
        
        return CommandExplanation(
            tool=tool,
            command=command,
            explanation=explanation,
            success=True
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error processing command: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    
    host = os.getenv("API_HOST", "localhost")
    port = int(os.getenv("API_PORT", 8001))  # Use different port
    
    uvicorn.run(app, host=host, port=port, reload=False)