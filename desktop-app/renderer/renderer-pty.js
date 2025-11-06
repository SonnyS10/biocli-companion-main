// BioCLI Companion - Professional Terminal with node-pty
console.log('🚀 Loading BioCLI Terminal with node-pty...');

let terminal;
let ptyProcess;
let commandHistory = [];
let currentCommand = '';

// Wait for DOM to be ready
document.addEventListener('DOMContentLoaded', () => {
    console.log('📄 DOM loaded, initializing professional terminal...');
    
    // Test basic functionality first
    testBasicElements();
    
    // Load the real terminal with node-pty
    loadProfessionalTerminal();
    
    // Setup UI interactions
    setupUI();
    
    console.log('✅ BioCLI Companion initialized!');
});

function testBasicElements() {
    console.log('🔍 Testing UI elements...');
    
    const terminalContainer = document.getElementById('terminal-container');
    const chatContainer = document.getElementById('chat-container');
    const userInput = document.getElementById('user-input');
    
    console.log('Terminal container:', terminalContainer ? '✅' : '❌');
    console.log('Chat container:', chatContainer ? '✅' : '❌');
    console.log('User input:', userInput ? '✅' : '❌');
    
    // Show loading state
    if (terminalContainer) {
        terminalContainer.innerHTML = '<div style="color: #7BC142; padding: 20px; text-align: center;">🧬 Initializing BioCLI Terminal...</div>';
    }
}

function loadProfessionalTerminal() {
    console.log('🖥️ Creating professional terminal with node-pty...');
    
    try {
        // Check if Terminal class is available
        if (typeof Terminal === 'undefined') {
            throw new Error('xterm.js Terminal not available');
        }
        
        console.log('✅ xterm.js loaded, creating terminal instance...');
        
        // Create terminal with professional settings
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
                white: '#ffffff',
                brightBlack: '#6c757d',
                brightRed: '#dc3545',
                brightGreen: '#28a745',
                brightYellow: '#ffc107',
                brightBlue: '#007bff',
                brightMagenta: '#6f42c1',
                brightCyan: '#17a2b8',
                brightWhite: '#ffffff'
            },
            fontFamily: '"Consolas", "Courier New", monospace',
            fontSize: 14,
            fontWeight: 'normal',
            fontWeightBold: 'bold',
            lineHeight: 1.2,
            letterSpacing: 0,
            rows: 30,
            cols: 100,
            scrollback: 1000,
            tabStopWidth: 8,
            bellSound: null,
            bellStyle: 'none'
        });

        // Mount terminal to container
        const terminalContainer = document.getElementById('terminal-container');
        if (!terminalContainer) {
            throw new Error('Terminal container not found');
        }
        
        terminal.open(terminalContainer);
        console.log('✅ Terminal mounted to DOM');
        
        // Show welcome message
        showWelcomeMessage();
        
        // Start the PTY process
        startPTYProcess();
        
    } catch (error) {
        console.error('❌ Failed to create terminal:', error);
        showTerminalError(error);
    }
}

function showWelcomeMessage() {
    console.log('🎨 Displaying welcome message...');
    
    terminal.write('\x1b[2J\x1b[H'); // Clear screen and move cursor to home
    terminal.write('\x1b[32m🧬 BioCLI Companion - Professional Terminal 🧬\x1b[0m\r\n');
    terminal.write('┌─────────────────────────────────────────────────────────┐\r\n');
    terminal.write('│ \x1b[36m🚀 Fully Interactive PowerShell Terminal\x1b[0m                 │\r\n');
    terminal.write('│ \x1b[33m• Execute real bioinformatics commands\x1b[0m                  │\r\n');
    terminal.write('│ \x1b[33m• Get instant AI explanations in sidebar\x1b[0m               │\r\n');
    terminal.write('│ \x1b[33m• Learn while you work with real tools!\x1b[0m                │\r\n');
    terminal.write('└─────────────────────────────────────────────────────────┘\r\n');
    terminal.write('\r\n\x1b[36m🔌 Connecting to PowerShell with node-pty...\x1b[0m\r\n\r\n');
}

function startPTYProcess() {
    console.log('🐚 Starting PTY process...');
    
    try {
        const pty = require('node-pty');
        
        // Create PTY process (Windows PowerShell)
        ptyProcess = pty.spawn('powershell.exe', [], {
            name: 'xterm-color',
            cols: terminal.cols,
            rows: terminal.rows,
            cwd: process.cwd(),
            env: process.env
        });
        
        console.log('✅ PTY process created successfully');
        terminal.write('\x1b[32m✅ Connected to PowerShell!\x1b[0m\r\n\r\n');
        
        // Handle PTY output
        ptyProcess.onData((data) => {
            terminal.write(data);
        });
        
        // Handle PTY exit
        ptyProcess.onExit((exitCode) => {
            console.log(`🚪 PTY process exited with code: ${exitCode.exitCode}`);
            terminal.write(`\r\n\r\n\x1b[33m🚪 Shell session ended (exit code: ${exitCode.exitCode})\x1b[0m\r\n`);
            terminal.write('\x1b[36m🔄 Restart the app to reconnect\x1b[0m\r\n');
        });
        
        // Handle terminal input -> PTY
        terminal.onData((data) => {
            if (ptyProcess) {
                ptyProcess.write(data);
                
                // Track commands for AI integration
                trackCommand(data);
            }
        });
        
        // Handle terminal resize
        terminal.onResize((size) => {
            if (ptyProcess) {
                ptyProcess.resize(size.cols, size.rows);
            }
        });
        
        console.log('✅ PTY communication established');
        
    } catch (error) {
        console.error('❌ Failed to start PTY process:', error);
        terminal.write(`\x1b[31m❌ Failed to start PowerShell: ${error.message}\x1b[0m\r\n`);
        terminal.write('\x1b[33m💡 Make sure node-pty is properly installed\x1b[0m\r\n');
        terminal.write('\x1b[36m📝 You can still use the AI chat on the right →\x1b[0m\r\n');
    }
}

function trackCommand(data) {
    // Build current command as user types
    if (data === '\r') {
        // Command entered (Enter key)
        if (currentCommand.trim()) {
            console.log('📝 Command entered:', currentCommand.trim());
            
            // Add to history
            commandHistory.push({
                command: currentCommand.trim(),
                timestamp: new Date()
            });
            
            // Check if it's a bioinformatics command
            if (isBioinformaticsCommand(currentCommand.trim())) {
                console.log('🧬 Bioinformatics command detected!');
                // TODO: Trigger AI explanation in sidebar
                addBioCommandNotification(currentCommand.trim());
            }
            
            currentCommand = '';
        }
    } else if (data === '\x7f' || data === '\b') {
        // Backspace
        currentCommand = currentCommand.slice(0, -1);
    } else if (data.length === 1 && data.charCodeAt(0) >= 32) {
        // Regular character
        currentCommand += data;
    }
}

function isBioinformaticsCommand(command) {
    const bioTools = [
        'bwa', 'samtools', 'fastqc', 'star', 'bcftools', 'bedtools',
        'bowtie2', 'tophat', 'cufflinks', 'gatk', 'picard', 'trimmomatic',
        'blast', 'blastn', 'blastp', 'blastx', 'tblastn', 'tblastx',
        'muscle', 'clustalw', 'mafft', 'hisat2', 'kallisto', 'salmon',
        'featurecounts', 'htseq-count', 'rsem', 'stringtie', 'cuffmerge'
    ];
    
    const firstWord = command.split(' ')[0].toLowerCase();
    return bioTools.includes(firstWord);
}

function addBioCommandNotification(command) {
    // Add visual notification in chat that we detected a bio command
    const chatContainer = document.getElementById('chat-container');
    if (chatContainer) {
        const notification = document.createElement('div');
        notification.className = 'bio-command-notification';
        notification.innerHTML = `
            <div style="background: rgba(123, 193, 66, 0.1); border-left: 3px solid #7BC142; padding: 10px; margin: 10px 0; border-radius: 5px;">
                <strong>🧬 Bioinformatics Command Detected!</strong><br>
                <code>${command}</code><br>
                <small>AI explanation coming soon...</small>
            </div>
        `;
        chatContainer.appendChild(notification);
        chatContainer.scrollTop = chatContainer.scrollHeight;
    }
}

function showTerminalError(error) {
    const terminalContainer = document.getElementById('terminal-container');
    if (terminalContainer) {
        terminalContainer.innerHTML = `
            <div style="color: #dc3545; padding: 20px; text-align: center;">
                <h3>⚠️ Terminal Initialization Error</h3>
                <p><strong>Error:</strong> ${error.message}</p>
                <p>Common solutions:</p>
                <ul style="text-align: left; display: inline-block;">
                    <li>Make sure node-pty is installed: <code>npm install node-pty</code></li>
                    <li>Check that you have Python and Visual Studio Build Tools</li>
                    <li>Restart the app after installation</li>
                </ul>
            </div>
        `;
    }
}

function setupUI() {
    console.log('🎛️ Setting up UI interactions...');
    
    // Chat input functionality
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
    
    // Control buttons
    setupControlButtons();
}

function handleChatMessage() {
    const userInput = document.getElementById('user-input');
    const message = userInput.value.trim();
    
    if (message) {
        console.log('💬 Chat message:', message);
        addChatMessage('user', message);
        
        // Demo AI response
        setTimeout(() => {
            addChatMessage('ai', `I can help explain bioinformatics commands! Try running commands like "bwa", "samtools", or "fastqc" in the terminal and I'll provide detailed explanations.`);
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
                const isHidden = aiPanel.style.display === 'none';
                aiPanel.style.display = isHidden ? 'flex' : 'none';
                
                // Resize terminal after sidebar toggle
                if (terminal && ptyProcess) {
                    setTimeout(() => {
                        const newCols = Math.floor(terminal.element.clientWidth / 9);
                        const newRows = Math.floor(terminal.element.clientHeight / 17);
                        terminal.resize(newCols, newRows);
                        ptyProcess.resize(newCols, newRows);
                    }, 100);
                }
            }
        });
    }
}

// Clean up on window close
window.addEventListener('beforeunload', () => {
    if (ptyProcess) {
        console.log('🧹 Cleaning up PTY process...');
        ptyProcess.kill();
    }
});

console.log('📝 Professional terminal script loaded');