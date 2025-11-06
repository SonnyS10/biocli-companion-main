// BioCLI Companion - Hybrid Terminal (Works Without node-pty)
console.log('🚀 Loading BioCLI Hybrid Terminal...');

let terminal;
let shellProcess;
let commandBuffer = '';
let commandHistory = [];

// Wait for DOM to be ready
document.addEventListener('DOMContentLoaded', () => {
    console.log('📄 DOM loaded, initializing hybrid terminal...');
    
    testBasicElements();
    loadHybridTerminal();
    setupUI();
    
    console.log('✅ BioCLI Companion initialized!');
});

function testBasicElements() {
    console.log('🔍 Testing UI elements...');
    
    const terminalContainer = document.getElementById('terminal-container');
    const chatContainer = document.getElementById('chat-container');
    
    console.log('Terminal container:', terminalContainer ? '✅' : '❌');
    console.log('Chat container:', chatContainer ? '✅' : '❌');
    
    if (terminalContainer) {
        terminalContainer.innerHTML = '<div style="color: #7BC142; padding: 20px; text-align: center;">🧬 Initializing BioCLI Hybrid Terminal...</div>';
    }
}

function loadHybridTerminal() {
    console.log('🖥️ Creating hybrid terminal...');
    
    try {
        if (typeof Terminal === 'undefined') {
            throw new Error('xterm.js Terminal not available');
        }
        
        // Create terminal with enhanced settings
        terminal = new Terminal({
            cursorBlink: true,
            theme: {
                background: '#000000',
                foreground: '#e0e0e0',
                cursor: '#7BC142',
                selection: '#2D5A87',
                black: '#000000',
                red: '#dc3545',
                green: '#28a745',
                yellow: '#ffc107',
                blue: '#007bff',
                magenta: '#6f42c1',
                cyan: '#17a2b8',
                white: '#ffffff'
            },
            fontFamily: '"Consolas", "Courier New", monospace',
            fontSize: 14,
            rows: 30,
            cols: 100,
            scrollback: 1000
        });

        const terminalContainer = document.getElementById('terminal-container');
        terminal.open(terminalContainer);
        console.log('✅ Terminal mounted to DOM');
        
        showWelcomeMessage();
        startHybridShell();
        
    } catch (error) {
        console.error('❌ Failed to create terminal:', error);
        showTerminalError(error);
    }
}

function showWelcomeMessage() {
    terminal.write('\x1b[2J\x1b[H'); // Clear screen
    terminal.write('\x1b[32m🧬 BioCLI Companion - Hybrid Terminal 🧬\x1b[0m\r\n');
    terminal.write('┌─────────────────────────────────────────────────────────┐\r\n');
    terminal.write('│ \x1b[36m🚀 Enhanced PowerShell Integration\x1b[0m                    │\r\n');
    terminal.write('│ \x1b[33m• Execute bioinformatics commands\x1b[0m                     │\r\n');
    terminal.write('│ \x1b[33m• Smart command detection & AI integration\x1b[0m            │\r\n');
    terminal.write('│ \x1b[33m• Works without node-pty dependencies!\x1b[0m               │\r\n');
    terminal.write('└─────────────────────────────────────────────────────────┘\r\n');
    terminal.write('\r\n\x1b[36m🔌 Connecting to PowerShell...\x1b[0m\r\n\r\n');
}

function startHybridShell() {
    console.log('🐚 Starting hybrid shell...');
    
    try {
        const { spawn } = require('child_process');
        
        // Start PowerShell with better configuration
        shellProcess = spawn('powershell.exe', [
            '-NoLogo',
            '-NoExit', 
            '-Command', 
            'Write-Host "BioCLI Ready!" -ForegroundColor Green; $Host.UI.RawUI.WindowTitle="BioCLI-PowerShell"'
        ], {
            cwd: process.cwd(),
            env: process.env,
            shell: false
        });
        
        console.log('✅ PowerShell started');
        terminal.write('\x1b[32m✅ Connected to Enhanced PowerShell!\x1b[0m\r\n\r\n');
        
        // Enhanced output handling
        shellProcess.stdout.on('data', (data) => {
            const output = data.toString();
            terminal.write(output);
            
            // Look for command completion patterns
            if (output.includes('PS ') || output.includes('> ')) {
                // Command likely completed, check for bio commands
                checkLastCommand();
            }
        });
        
        shellProcess.stderr.on('data', (data) => {
            const error = data.toString();
            terminal.write(`\x1b[31m${error}\x1b[0m`);
        });
        
        shellProcess.on('exit', (code) => {
            console.log(`🚪 Shell exited with code: ${code}`);
            terminal.write(`\r\n\x1b[33m🚪 Shell session ended\x1b[0m\r\n`);
        });
        
        // Enhanced input handling
        setupEnhancedInput();
        
    } catch (error) {
        console.error('❌ Failed to start shell:', error);
        terminal.write(`\x1b[31m❌ Failed to start PowerShell: ${error.message}\x1b[0m\r\n`);
    }
}

function setupEnhancedInput() {
    let currentLine = '';
    
    terminal.onData((data) => {
        if (!shellProcess || shellProcess.killed) return;
        
        const code = data.charCodeAt(0);
        
        if (data === '\r') {
            // Enter pressed - command submitted
            if (currentLine.trim()) {
                const command = currentLine.trim();
                console.log('📝 Command entered:', command);
                
                // Add to history
                commandHistory.push({
                    command: command,
                    timestamp: new Date()
                });
                
                // Check for bioinformatics commands
                if (isBioinformaticsCommand(command)) {
                    console.log('🧬 Bioinformatics command detected!');
                    addBioCommandNotification(command);
                }
                
                currentLine = '';
            }
            
            // Send to shell
            shellProcess.stdin.write(data);
            
        } else if (data === '\x7f' || data === '\b') {
            // Backspace
            if (currentLine.length > 0) {
                currentLine = currentLine.slice(0, -1);
            }
            shellProcess.stdin.write(data);
            
        } else if (code >= 32 || data === '\t') {
            // Regular character or tab
            currentLine += data;
            shellProcess.stdin.write(data);
            
        } else {
            // Other control characters
            shellProcess.stdin.write(data);
        }
    });
}

function checkLastCommand() {
    if (commandHistory.length > 0) {
        const lastCommand = commandHistory[commandHistory.length - 1];
        const timeSinceCommand = Date.now() - lastCommand.timestamp.getTime();
        
        // If command was recent and bioinformatics-related
        if (timeSinceCommand < 2000 && isBioinformaticsCommand(lastCommand.command)) {
            // Add enhanced notification
            addEnhancedBioNotification(lastCommand.command);
        }
    }
}

function isBioinformaticsCommand(command) {
    const bioTools = [
        'bwa', 'samtools', 'fastqc', 'star', 'bcftools', 'bedtools',
        'bowtie2', 'tophat', 'cufflinks', 'gatk', 'picard', 'trimmomatic',
        'blast', 'blastn', 'blastp', 'blastx', 'tblastn', 'tblastx',
        'muscle', 'clustalw', 'mafft', 'hisat2', 'kallisto', 'salmon',
        'featurecounts', 'htseq-count', 'rsem', 'stringtie'
    ];
    
    const firstWord = command.split(' ')[0].toLowerCase();
    return bioTools.includes(firstWord);
}

function addBioCommandNotification(command) {
    const chatContainer = document.getElementById('chat-container');
    if (chatContainer) {
        const notification = document.createElement('div');
        notification.className = 'bio-command-notification';
        notification.innerHTML = `
            <div style="background: rgba(123, 193, 66, 0.1); border-left: 3px solid #7BC142; padding: 12px; margin: 8px 0; border-radius: 6px;">
                <div style="font-weight: bold; color: #7BC142; margin-bottom: 6px;">
                    🧬 Bioinformatics Command Detected!
                </div>
                <code style="background: rgba(0,0,0,0.3); padding: 4px 8px; border-radius: 4px;">${command}</code>
                <div style="margin-top: 8px; font-size: 12px; color: #b0b0b0;">
                    💡 AI explanation will be available once backend is connected
                </div>
            </div>
        `;
        chatContainer.appendChild(notification);
        chatContainer.scrollTop = chatContainer.scrollHeight;
    }
}

function addEnhancedBioNotification(command) {
    // More detailed notification for completed commands
    addChatMessage('ai', `🧬 I detected you just ran: <code>${command}</code><br><br>This is a bioinformatics tool! Once we connect the AI backend, I'll provide detailed explanations of what this command does, its parameters, and best practices.`);
}

function showTerminalError(error) {
    const terminalContainer = document.getElementById('terminal-container');
    if (terminalContainer) {
        terminalContainer.innerHTML = `
            <div style="color: #dc3545; padding: 20px;">
                <h3>⚠️ Terminal Error</h3>
                <p>${error.message}</p>
                <p>Using fallback terminal mode...</p>
            </div>
        `;
    }
}

// UI Functions (same as before)
function setupUI() {
    const userInput = document.getElementById('user-input');
    const sendButton = document.getElementById('send-button');
    
    if (userInput && sendButton) {
        sendButton.addEventListener('click', handleChatMessage);
        userInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                handleChatMessage();
            }
        });
    }
    
    setupControlButtons();
}

function handleChatMessage() {
    const userInput = document.getElementById('user-input');
    const message = userInput.value.trim();
    
    if (message) {
        addChatMessage('user', message);
        
        setTimeout(() => {
            if (isBioinformaticsCommand(message)) {
                addChatMessage('ai', `That's a bioinformatics command! Try running <code>${message}</code> in the terminal on the left, and I'll detect it automatically.`);
            } else {
                addChatMessage('ai', 'I can help explain bioinformatics commands! Try running tools like "bwa", "samtools", or "fastqc" in the terminal.');
            }
        }, 500);
        
        userInput.value = '';
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

function setupControlButtons() {
    // Clear chat
    const clearChatBtn = document.getElementById('clear-chat');
    if (clearChatBtn) {
        clearChatBtn.addEventListener('click', () => {
            const chatContainer = document.getElementById('chat-container');
            if (chatContainer) {
                const messages = chatContainer.querySelectorAll('.message, .bio-command-notification');
                messages.forEach(msg => msg.remove());
            }
        });
    }
    
    // Clear terminal
    const clearTerminalBtn = document.getElementById('clear-terminal');
    if (clearTerminalBtn) {
        clearTerminalBtn.addEventListener('click', () => {
            if (terminal) {
                terminal.clear();
                showWelcomeMessage();
            }
        });
    }
    
    // Toggle sidebar
    const toggleBtn = document.getElementById('toggle-sidebar');
    if (toggleBtn) {
        toggleBtn.addEventListener('click', () => {
            const aiPanel = document.getElementById('ai-panel');
            if (aiPanel) {
                aiPanel.style.display = aiPanel.style.display === 'none' ? 'flex' : 'none';
            }
        });
    }
}

// Clean up
window.addEventListener('beforeunload', () => {
    if (shellProcess && !shellProcess.killed) {
        shellProcess.kill();
    }
});

console.log('📝 Hybrid terminal script loaded');