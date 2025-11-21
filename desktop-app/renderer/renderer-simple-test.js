// BioCLI Companion - AI-Powered Terminal with Real Shell
console.log('🚀 Loading BioCLI Terminal with Shell Integration...');

const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');

let terminal;
let shellProcess = null;
let currentCommand = '';
let currentPath = process.env.USERPROFILE || process.env.HOME || 'C:\\';
let commandHistory = [];
let historyIndex = -1;

document.addEventListener('DOMContentLoaded', () => {
    console.log('📄 DOM loaded, creating terminal with shell...');
    
    initializeTerminal();
    setupBasicUI();
    
    console.log('✅ Terminal with shell ready!');
});

function initializeTerminal() {
    console.log('🖥️ Creating terminal with shell integration...');
    
    try {
        if (typeof Terminal === 'undefined') {
            throw new Error('xterm.js not loaded');
        }
        
        // Create terminal
        terminal = new Terminal({
            cursorBlink: true,
            theme: {
                background: '#000000',
                foreground: '#e0e0e0',
                cursor: '#7BC142'
            },
            fontFamily: 'Consolas, monospace',
            fontSize: 14,
            rows: 25,
            cols: 80
        });

        const container = document.getElementById('terminal-container');
        if (!container) {
            throw new Error('Terminal container not found');
        }
        
        terminal.open(container);
        console.log('✅ Terminal created and mounted');
        
        // Show welcome message
        showWelcomeMessage();
        
        // Start PowerShell
        startShell();
        
        // Setup input handling
        setupTerminalInput();
        
    } catch (error) {
        console.error('❌ Terminal creation failed:', error);
        showError(error);
    }
}

function showWelcomeMessage() {
    terminal.write('🧬 BioCLI Companion - AI-Powered Terminal 🧬\r\n');
    terminal.write('═══════════════════════════════════════\r\n');
    terminal.write('🚀 Real PowerShell integration active!\r\n');
    terminal.write('💡 Type commands like: dir, echo, cd, ls\r\n');
    terminal.write('🧬 Bioinformatics tools: bwa, samtools, fastqc\r\n');
    terminal.write('🤖 AI Backend connected!\r\n');
    terminal.write('📝 Type "demo" to see AI formatting\r\n\r\n');
}

function startShell() {
    console.log('🐚 Starting PowerShell...');
    
    try {
        // Get initial directory
        currentPath = process.env.USERPROFILE || process.env.HOME || 'C:\\';
        
        // Spawn PowerShell process
        shellProcess = spawn('powershell.exe', [
            '-NoLogo', 
            '-NoProfile',
            '-NoExit',
            '-Command', '-'
        ], {
            cwd: currentPath,
            env: process.env,
            shell: false
        });
        
        let outputBuffer = '';
        let isFirstOutput = true;
        
        // Handle shell output
        shellProcess.stdout.on('data', (data) => {
            const output = data.toString();
            
            // Skip initial PowerShell startup messages
            if (isFirstOutput) {
                isFirstOutput = false;
                return;
            }
            
            // Check if output contains path (PowerShell prompt)
            const pathMatch = output.match(/PS ([A-Z]:[^\>]+)>/);
            if (pathMatch) {
                currentPath = pathMatch[1];
            }
            
            // Write output to terminal
            terminal.write(output.replace(/\n/g, '\r\n'));
        });
        
        // Handle errors (including "command not found")
        shellProcess.stderr.on('data', (data) => {
            const error = data.toString();
            // Show errors in red
            terminal.write(`\x1b[31m${error.replace(/\n/g, '\r\n')}\x1b[0m`);
        });
        
        shellProcess.on('error', (error) => {
            console.error('❌ Shell process error:', error);
            terminal.write(`\r\n\x1b[31mShell error: ${error.message}\x1b[0m\r\n`);
        });
        
        shellProcess.on('exit', (code) => {
            console.log(`🚪 Shell exited with code ${code}`);
            terminal.write(`\r\n\x1b[33mShell exited with code ${code}. Press F5 to restart.\x1b[0m\r\n`);
            shellProcess = null;
        });
        
        console.log('✅ PowerShell started successfully');
        
        // Give shell a moment to initialize, then show prompt
        setTimeout(() => {
            showPrompt();
        }, 500);
        
    } catch (error) {
        console.error('❌ Failed to start shell:', error);
        terminal.write(`\r\n\x1b[31mFailed to start PowerShell: ${error.message}\x1b[0m\r\n`);
        terminal.write('Falling back to built-in commands mode.\r\n');
        shellProcess = null;
        showPrompt();
    }
}

function setupTerminalInput() {
    terminal.onData((data) => {
        const char = data;
        
        if (char === '\r') {
            // Enter pressed
            terminal.write('\r\n');
            
            if (currentCommand.trim()) {
                // Add to history
                commandHistory.push(currentCommand.trim());
                historyIndex = commandHistory.length;
                
                handleTerminalCommand(currentCommand.trim());
            } else {
                // Empty command, just show new prompt
                showPrompt();
            }
            
            currentCommand = '';
            
        } else if (char === '\t') {
            // Tab pressed - autocomplete
            handleTabComplete();
            
        } else if (char === '\x7f' || char === '\b') {
            // Backspace
            if (currentCommand.length > 0) {
                currentCommand = currentCommand.slice(0, -1);
                terminal.write('\b \b');
            }
            
        } else if (char === '\x03') {
            // Ctrl+C
            terminal.write('^C\r\n');
            currentCommand = '';
            showPrompt();
            
        } else if (char === '\x1b[A') {
            // Up arrow - previous command
            handleHistoryUp();
            
        } else if (char === '\x1b[B') {
            // Down arrow - next command
            handleHistoryDown();
            
        } else if (char >= ' ' && char.length === 1) {
            // Regular character
            currentCommand += char;
            terminal.write(char);
        }
    });
}

function handleTabComplete() {
    if (!currentCommand.trim()) return;
    
    const parts = currentCommand.split(' ');
    const lastPart = parts[parts.length - 1];
    
    // Try to autocomplete file/directory names
    try {
        const searchPath = lastPart.includes('\\') || lastPart.includes('/') 
            ? path.dirname(lastPart) 
            : currentPath;
        
        const searchTerm = lastPart.includes('\\') || lastPart.includes('/') 
            ? path.basename(lastPart)
            : lastPart;
        
        // Get directory contents
        const files = fs.readdirSync(searchPath);
        
        // Find matches
        const matches = files.filter(file => 
            file.toLowerCase().startsWith(searchTerm.toLowerCase())
        );
        
        if (matches.length === 1) {
            // Single match - autocomplete it
            const completion = matches[0];
            const toAdd = completion.substring(searchTerm.length);
            
            currentCommand += toAdd;
            terminal.write(toAdd);
            
            // If it's a directory, add a backslash
            const fullPath = path.join(searchPath, completion);
            try {
                if (fs.statSync(fullPath).isDirectory()) {
                    currentCommand += '\\';
                    terminal.write('\\');
                }
            } catch (e) {
                // Ignore stat errors
            }
            
        } else if (matches.length > 1) {
            // Multiple matches - show them
            terminal.write('\r\n');
            matches.forEach(match => {
                terminal.write(`  ${match}\r\n`);
            });
            showPrompt();
            terminal.write(currentCommand);
        }
        
    } catch (error) {
        // Silently fail autocomplete
        console.log('Autocomplete error:', error.message);
    }
}

function handleHistoryUp() {
    if (historyIndex > 0) {
        historyIndex--;
        replaceCurrentCommand(commandHistory[historyIndex]);
    }
}

function handleHistoryDown() {
    if (historyIndex < commandHistory.length - 1) {
        historyIndex++;
        replaceCurrentCommand(commandHistory[historyIndex]);
    } else if (historyIndex === commandHistory.length - 1) {
        historyIndex++;
        replaceCurrentCommand('');
    }
}

function replaceCurrentCommand(newCommand) {
    // Clear current command
    for (let i = 0; i < currentCommand.length; i++) {
        terminal.write('\b \b');
    }
    
    // Write new command
    currentCommand = newCommand;
    terminal.write(newCommand);
}

function handleTerminalCommand(command) {
    console.log('📝 Command entered:', command);
    
    // Check for built-in demo/help commands first
    if (command.toLowerCase() === 'help') {
        showHelp();
        return;
    }
    
    if (command.toLowerCase() === 'demo') {
        terminal.write('🎨 Showing demo AI response...\r\n');
        showDemoResponse();
        showPrompt();
        return;
    }
    
    if (command.toLowerCase() === 'clear') {
        terminal.clear();
        showWelcomeMessage();
        showPrompt();
        return;
    }
    
    // Translate Unix-style commands to PowerShell equivalents
    command = translateUnixCommand(command);
    
    // Track cd commands to update current path
    if (command.toLowerCase().startsWith('cd ')) {
        const targetPath = command.substring(3).trim();
        try {
            const newPath = path.resolve(currentPath, targetPath);
            if (fs.existsSync(newPath) && fs.statSync(newPath).isDirectory()) {
                currentPath = newPath;
            }
        } catch (e) {
            // Let PowerShell handle the error
        }
    }
    
    // Check for bioinformatics commands
    const bioTools = ['bwa', 'samtools', 'fastqc', 'star', 'bcftools', 'bowtie2', 'gatk', 'picard', 'bedtools'];
    const commandLower = command.toLowerCase().split(' ')[0];
    
    if (bioTools.includes(commandLower)) {
        terminal.write(`🧬 Bioinformatics tool detected: ${commandLower}\r\n`);
        getAIExplanation(command);
    }
    
    // Execute in shell if available
    if (shellProcess && !shellProcess.killed) {
        try {
            // Write command to PowerShell
            shellProcess.stdin.write(command + '\n');
            
            // Show prompt after command completes (with delay for output)
            setTimeout(() => {
                showPrompt();
            }, 100);
            
        } catch (error) {
            console.error('❌ Error executing command:', error);
            terminal.write(`\x1b[31mError: ${error.message}\x1b[0m\r\n`);
            showPrompt();
        }
    } else {
        // Fallback: show that shell is not available
        terminal.write(`\x1b[33mPowerShell not available.\x1b[0m\r\n`);
        terminal.write(`Command entered: ${command}\r\n`);
        showPrompt();
    }
}

// Translate common Unix commands to PowerShell equivalents
function translateUnixCommand(command) {
    const parts = command.split(' ');
    const cmd = parts[0].toLowerCase();
    const args = parts.slice(1).join(' ');
    
    const translations = {
        // ls variations
        'ls -la': 'Get-ChildItem -Force',
        'ls -al': 'Get-ChildItem -Force',
        'ls -a': 'Get-ChildItem -Force',
        'ls -l': 'Get-ChildItem',
        
        // grep
        'grep': 'Select-String',
        
        // cat
        'cat': 'Get-Content',
        
        // touch
        'touch': 'New-Item -ItemType File',
        
        // rm variations
        'rm -rf': 'Remove-Item -Recurse -Force',
        'rm -r': 'Remove-Item -Recurse',
        'rm -f': 'Remove-Item -Force',
        
        // cp variations
        'cp -r': 'Copy-Item -Recurse',
        
        // mv
        'mv': 'Move-Item',
        
        // mkdir
        'mkdir -p': 'New-Item -ItemType Directory -Force',
        
        // pwd
        'pwd': 'Get-Location',
        
        // which
        'which': 'Get-Command',
        
        // find
        'find': 'Get-ChildItem -Recurse',
        
        // head
        'head': 'Get-Content -TotalCount 10',
        
        // tail
        'tail': 'Get-Content -Tail 10',
        
        // wc
        'wc -l': 'Measure-Object -Line'
    };
    
    // Check for exact matches first
    if (translations[command.toLowerCase()]) {
        const translated = translations[command.toLowerCase()];
        terminal.write(`\x1b[90m(Translated: ${command} → ${translated})\x1b[0m\r\n`);
        return translated;
    }
    
    // Check for command with args
    if (translations[cmd]) {
        const translated = `${translations[cmd]} ${args}`;
        terminal.write(`\x1b[90m(Translated: ${command} → ${translated})\x1b[0m\r\n`);
        return translated;
    }
    
    // No translation needed
    return command;
}

function showPrompt() {
    if (shellProcess && !shellProcess.killed) {
        // Show PowerShell-style prompt with path
        terminal.write(`\x1b[36mPS ${currentPath}>\x1b[0m `);
    } else {
        terminal.write('biocli> ');
    }
}

function showHelp() {
    terminal.write('\r\n🧬 BioCLI Companion - Command Help\r\n');
    terminal.write('═══════════════════════════════════════\r\n');
    terminal.write('Built-in Commands:\r\n');
    terminal.write('  help  - Show this help\r\n');
    terminal.write('  demo  - Show demo AI response\r\n');
    terminal.write('  clear - Clear terminal\r\n');
    terminal.write('\r\nPowerShell Commands:\r\n');
    terminal.write('  Any PowerShell command (dir, ls, cd, Get-Date, etc.)\r\n');
    terminal.write('  Current path shown in prompt!\r\n');
    terminal.write('  Errors will display in red\r\n');
    terminal.write('\r\nBioinformatics Tools (AI-enhanced):\r\n');
    terminal.write('  bwa, samtools, fastqc, star, bcftools, bedtools\r\n');
    terminal.write('  picard, bowtie2, gatk\r\n');
    terminal.write('  → These will trigger AI explanations!\r\n');
    terminal.write('\r\nKeyboard Shortcuts:\r\n');
    terminal.write('  Tab         - Autocomplete files/folders\r\n');
    terminal.write('  Up/Down     - Command history\r\n');
    terminal.write('  Ctrl+C      - Cancel current input\r\n');
    terminal.write('  Backspace   - Delete character\r\n\r\n');
    
    showPrompt();
}

// NEW: AI Backend Integration
const API_BASE_URL = 'http://127.0.0.1:8000';

// Demo response to show formatting without API
function showDemoResponse() {
    const demoData = {
        tool: 'BWA',
        command: 'bwa mem',
        explanation: `## Overview
**BWA MEM** is a fast and accurate alignment algorithm for mapping sequencing reads to a reference genome. It's part of the Burrows-Wheeler Aligner (BWA) suite.

## Key Features
- **Ultra-fast alignment** for short-read sequencing data
- Supports paired-end and single-end reads
- Works with Illumina, 454, and Ion Torrent data
- *Ideal for reads between 70bp - 1Mbp*

## Basic Usage

\`\`\`bash
bwa mem reference.fa reads.fq > output.sam
\`\`\`

## Common Parameters
- \`-t\` - Number of threads (default: 1)
- \`-M\` - Mark shorter split hits as secondary
- \`-R\` - Add read group information

## Example: Paired-End Alignment

\`\`\`bash
bwa mem -t 8 -M reference.fa read1.fq read2.fq > aligned.sam
\`\`\`

## Best Practices
1. **Always index your reference** first with \`bwa index\`
2. Use multiple threads (\`-t\`) for faster processing
3. Add read group info for downstream analysis
4. Convert SAM to BAM to save space

**Note:** BWA MEM is the recommended algorithm for reads longer than 70bp.`
    };
    
    addAIExplanation('bwa mem', demoData);
    terminal.write(`✅ Demo response displayed! Check the chat panel →\r\n`);
}

async function getAIExplanation(command) {
    const chatContainer = document.getElementById('chat-container');
    
    // Add loading message
    const loadingDiv = document.createElement('div');
    loadingDiv.id = 'ai-loading';
    loadingDiv.innerHTML = `
        <div style="background: rgba(123, 193, 66, 0.1); border-left: 3px solid #7BC142; padding: 12px; margin: 10px 0; border-radius: 6px;">
            <div style="font-weight: bold; color: #7BC142; margin-bottom: 8px;">
                🤖 AI is thinking...
            </div>
            <code style="background: rgba(0,0,0,0.3); padding: 4px 8px; border-radius: 4px; font-family: Consolas;">${command}</code>
            <div style="margin-top: 8px; font-size: 12px; color: #b0b0b0;">
                ⏳ Generating explanation...
            </div>
        </div>
    `;
    chatContainer.appendChild(loadingDiv);
    chatContainer.scrollTop = chatContainer.scrollHeight;
    
    try {
        console.log('📡 Calling AI backend with command:', command);
        
        // Call your FastAPI backend
        const response = await fetch(`${API_BASE_URL}/api/explain`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ command: command })
        });
        
        const data = await response.json();
        console.log('📦 Backend response:', data);
        
        // Remove loading message
        loadingDiv.remove();
        
        if (!response.ok || !data.success) {
            // Backend returned an error
            const errorMessage = data.error_message || data.explanation || 'Unknown error';
            throw new Error(errorMessage);
        }
        
        console.log('✅ AI response received:', data);
        
        // Add AI explanation
        addAIExplanation(command, data);
        
        // Also show in terminal
        terminal.write(`✅ AI explanation received! Check the chat panel →\r\n`);
        
    } catch (error) {
        console.error('❌ AI request failed:', error);
        
        // Remove loading message
        loadingDiv.remove();
        
        // Show detailed error
        const errorDiv = document.createElement('div');
        
        // Check for specific error types
        const isQuotaError = error.message.includes('quota') || error.message.includes('insufficient_quota');
        const isRateLimitError = error.message.includes('rate limit') || error.message.includes('429');
        
        let troubleshootingHTML = '';
        
        if (isQuotaError || isRateLimitError) {
            troubleshootingHTML = `
                <div style="margin-top: 12px; font-size: 12px; color: #b0b0b0; border-top: 1px solid rgba(255,255,255,0.1); padding-top: 8px;">
                    <strong style="color: #ffc107;">⚠️ OpenAI Quota Exceeded</strong><br>
                    <div style="margin-top: 8px;">
                        Your OpenAI API has exceeded its quota. To fix this:<br>
                        1. Check your OpenAI billing at <a href="https://platform.openai.com/account/billing" target="_blank" style="color: #7BC142;">platform.openai.com/account/billing</a><br>
                        2. Add payment method or upgrade plan<br>
                        3. Or use a different API key in backend/.env
                    </div>
                    <div style="margin-top: 12px; padding: 8px; background: rgba(123, 193, 66, 0.1); border-radius: 4px;">
                        <strong style="color: #7BC142;">✅ Good News!</strong><br>
                        Your app is working perfectly! The connection between your desktop app and backend is successful. You just need to add OpenAI credits.
                    </div>
                </div>
            `;
        } else {
            troubleshootingHTML = `
                <div style="margin-top: 12px; font-size: 12px; color: #b0b0b0; border-top: 1px solid rgba(255,255,255,0.1); padding-top: 8px;">
                    <strong>Troubleshooting:</strong><br>
                    • Check if backend is running: <code style="background: rgba(0,0,0,0.3); padding: 2px 6px; border-radius: 3px;">uvicorn main:app --reload</code><br>
                    • Verify OpenAI API key is set in backend/.env<br>
                    • Check backend terminal for error messages<br>
                    • Backend URL: ${API_BASE_URL}
                </div>
            `;
        }
        
        errorDiv.innerHTML = `
            <div style="background: rgba(220, 53, 69, 0.1); border-left: 3px solid #dc3545; padding: 12px; margin: 10px 0; border-radius: 6px;">
                <div style="font-weight: bold; color: #dc3545; margin-bottom: 8px;">
                    ❌ AI Request Failed
                </div>
                <div style="font-size: 13px; color: #e0e0e0; margin-bottom: 8px;">
                    ${error.message}
                </div>
                ${troubleshootingHTML}
            </div>
        `;
        chatContainer.appendChild(errorDiv);
        chatContainer.scrollTop = chatContainer.scrollHeight;
        
        terminal.write(`❌ AI request failed: ${error.message}\r\n`);
    }
}

// Simple markdown-like formatter for AI responses
function formatMarkdown(text) {
    if (!text) return 'No explanation available';
    
    let formatted = text;
    
    // Code blocks (```code```)
    formatted = formatted.replace(/```(\w+)?\n([\s\S]*?)```/g, (match, lang, code) => {
        return `<pre style="background: rgba(0,0,0,0.5); padding: 12px; border-radius: 6px; overflow-x: auto; border-left: 3px solid #7BC142; margin: 10px 0;"><code style="color: #7BC142; font-family: Consolas, monospace; font-size: 13px;">${escapeHtml(code.trim())}</code></pre>`;
    });
    
    // Inline code (`code`)
    formatted = formatted.replace(/`([^`]+)`/g, '<code style="background: rgba(123, 193, 66, 0.2); padding: 2px 6px; border-radius: 3px; color: #7BC142; font-family: Consolas, monospace; font-size: 13px;">$1</code>');
    
    // Bold (**text**)
    formatted = formatted.replace(/\*\*([^*]+)\*\*/g, '<strong style="color: #7BC142; font-weight: 600;">$1</strong>');
    
    // Italic (*text*)
    formatted = formatted.replace(/\*([^*]+)\*/g, '<em style="color: #b0b0b0; font-style: italic;">$1</em>');
    
    // Headers (### Header)
    formatted = formatted.replace(/^### (.+)$/gm, '<h3 style="color: #7BC142; font-size: 16px; font-weight: 600; margin: 16px 0 8px 0; border-bottom: 1px solid rgba(123, 193, 66, 0.3); padding-bottom: 4px;">$1</h3>');
    formatted = formatted.replace(/^## (.+)$/gm, '<h2 style="color: #7BC142; font-size: 18px; font-weight: 600; margin: 18px 0 10px 0; border-bottom: 2px solid rgba(123, 193, 66, 0.4); padding-bottom: 6px;">$1</h2>');
    
    // Bullet points (- item or * item)
    formatted = formatted.replace(/^[•\-\*] (.+)$/gm, '<div style="margin-left: 20px; margin-bottom: 6px; color: #e0e0e0;"><span style="color: #7BC142; margin-right: 8px;">•</span>$1</div>');
    
    // Numbered lists (1. item)
    formatted = formatted.replace(/^(\d+)\. (.+)$/gm, '<div style="margin-left: 20px; margin-bottom: 6px; color: #e0e0e0;"><span style="color: #7BC142; margin-right: 8px; font-weight: 600;">$1.</span>$2</div>');
    
    // Line breaks
    formatted = formatted.replace(/\n\n/g, '<br><br>');
    formatted = formatted.replace(/\n/g, '<br>');
    
    return formatted;
}

function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

function addAIExplanation(command, data) {
    const chatContainer = document.getElementById('chat-container');
    
    const explanationDiv = document.createElement('div');
    explanationDiv.className = 'message ai';
    
    const time = new Date().toLocaleTimeString();
    
    // Format the explanation with markdown-like styling
    const formattedExplanation = formatMarkdown(data.explanation);
    
    explanationDiv.innerHTML = `
        <div class="message-header">
            <span class="message-sender ai">🤖 AI Assistant</span>
            <span class="message-time">${time}</span>
        </div>
        <div class="message-content">
            <!-- Command Header -->
            <div style="background: linear-gradient(135deg, rgba(123, 193, 66, 0.15), rgba(123, 193, 66, 0.05)); border-left: 4px solid #7BC142; padding: 14px; border-radius: 8px; margin-bottom: 16px; box-shadow: 0 2px 8px rgba(0,0,0,0.2);">
                <div style="display: flex; align-items: center; gap: 12px; margin-bottom: 10px;">
                    <span style="font-size: 24px;">🧬</span>
                    <div>
                        <div style="font-size: 11px; color: #888; text-transform: uppercase; letter-spacing: 1px; margin-bottom: 4px;">Command Detected</div>
                        <code style="background: rgba(0,0,0,0.4); padding: 4px 10px; border-radius: 4px; font-family: Consolas, monospace; font-size: 14px; color: #7BC142; font-weight: 500;">${command}</code>
                    </div>
                </div>
                <div style="display: flex; align-items: center; gap: 8px; padding-top: 8px; border-top: 1px solid rgba(123, 193, 66, 0.2);">
                    <span style="font-size: 18px;">🔧</span>
                    <div>
                        <span style="font-size: 11px; color: #888; text-transform: uppercase; letter-spacing: 1px;">Tool:</span>
                        <span style="color: #7BC142; font-weight: 600; margin-left: 6px; font-size: 14px;">${data.tool || 'Unknown'}</span>
                    </div>
                </div>
            </div>
            
            <!-- AI Explanation Content -->
            <div style="
                line-height: 1.8; 
                color: #e0e0e0; 
                font-size: 14px;
                background: rgba(0,0,0,0.2);
                padding: 16px;
                border-radius: 8px;
                border: 1px solid rgba(123, 193, 66, 0.1);
            ">
                ${formattedExplanation}
            </div>
            
            <!-- Footer -->
            <div style="margin-top: 12px; padding-top: 12px; border-top: 1px solid rgba(255,255,255,0.1); font-size: 11px; color: #666; text-align: right;">
                <span style="margin-right: 8px;">⚡ Powered by AI</span>
                <span>${time}</span>
            </div>
        </div>
    `;
    
    chatContainer.appendChild(explanationDiv);
    chatContainer.scrollTop = chatContainer.scrollHeight;
}

function addBioCommandToChat(command) {
    const chatContainer = document.getElementById('chat-container');
    if (chatContainer) {
        const notification = document.createElement('div');
        notification.innerHTML = `
            <div style="background: rgba(123, 193, 66, 0.1); border-left: 3px solid #7BC142; padding: 12px; margin: 10px 0; border-radius: 6px;">
                <div style="font-weight: bold; color: #7BC142; margin-bottom: 8px;">
                    🧬 Bioinformatics Command Detected!
                </div>
                <code style="background: rgba(0,0,0,0.3); padding: 4px 8px; border-radius: 4px; font-family: Consolas;">${command}</code>
                <div style="margin-top: 8px; font-size: 12px; color: #b0b0b0;">
                    ✨ AI explanation will appear here once backend is connected
                </div>
            </div>
        `;
        chatContainer.appendChild(notification);
        chatContainer.scrollTop = chatContainer.scrollHeight;
    }
}

function showError(error) {
    const container = document.getElementById('terminal-container');
    if (container) {
        container.innerHTML = `
            <div style="color: #dc3545; padding: 20px;">
                <h3>⚠️ Terminal Error</h3>
                <p>Error: ${error.message}</p>
                <p>Check browser console for details.</p>
            </div>
        `;
    }
}

function setupBasicUI() {
    console.log('🎛️ Setting up basic UI...');
    
    // Chat input
    const userInput = document.getElementById('user-input');
    const sendButton = document.getElementById('send-button');
    
    if (userInput && sendButton) {
        sendButton.addEventListener('click', () => {
            const message = userInput.value.trim();
            if (message) {
                addChatMessage('user', message);
                addChatMessage('ai', 'Terminal UI is working! Try typing commands in the terminal on the left. I can detect bioinformatics commands like "bwa", "samtools", etc.');
                userInput.value = '';
            }
        });
        
        userInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                sendButton.click();
            }
        });
    }
    
    // Clear chat button
    const clearChat = document.getElementById('clear-chat');
    if (clearChat) {
        clearChat.addEventListener('click', () => {
            const chatContainer = document.getElementById('chat-container');
            if (chatContainer) {
                const messages = chatContainer.querySelectorAll('.message');
                messages.forEach(msg => msg.remove());
            }
        });
    }
    
    // Clear terminal button
    const clearTerminal = document.getElementById('clear-terminal');
    if (clearTerminal) {
        clearTerminal.addEventListener('click', () => {
            if (terminal) {
                showSimplePrompt();
            }
        });
    }
}

function addChatMessage(sender, content) {
    const chatContainer = document.getElementById('chat-container');
    if (!chatContainer) return;
    
    const messageDiv = document.createElement('div');
    messageDiv.className = `message ${sender}`;
    
    const time = new Date().toLocaleTimeString();
    
    messageDiv.innerHTML = `
        <div class="message-header">
            <span class="message-sender ${sender}">${sender === 'user' ? 'You' : 'AI Assistant'}</span>
            <span class="message-time">${time}</span>
        </div>
        <div class="message-content">${content}</div>
    `;
    
    chatContainer.appendChild(messageDiv);
    chatContainer.scrollTop = chatContainer.scrollHeight;
}

console.log('📝 Simple terminal script loaded');