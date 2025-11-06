// BioCLI Companion - Simple Terminal Integration
console.log('🚀 Renderer starting...');

// Wait for DOM to be ready
document.addEventListener('DOMContentLoaded', () => {
    console.log('📄 DOM loaded, initializing app...');
    
    // Test basic functionality first
    testBasicElements();
    
    // Try to load xterm
    loadTerminal();
    
    // Setup basic UI interactions
    setupBasicUI();
});

function testBasicElements() {
    console.log('🔍 Testing basic elements...');
    
    // Test if our HTML elements exist
    const terminalContainer = document.getElementById('terminal-container');
    const chatContainer = document.getElementById('chat-container');
    const userInput = document.getElementById('user-input');
    
    console.log('Terminal container:', terminalContainer ? '✅ Found' : '❌ Missing');
    console.log('Chat container:', chatContainer ? '✅ Found' : '❌ Missing');
    console.log('User input:', userInput ? '✅ Found' : '❌ Missing');
    
    // Update status
    if (terminalContainer) {
        terminalContainer.innerHTML = '<div style="color: #7BC142; padding: 20px;">🧬 Terminal loading...</div>';
    }
}

let terminal;
let shellProcess;

function loadTerminal() {
    console.log('🖥️ Attempting to load interactive terminal...');
    
    try {
        // Check if Terminal class is available
        if (typeof Terminal !== 'undefined') {
            console.log('✅ Terminal class found, creating interactive terminal...');
            
            // Create terminal instance
            terminal = new Terminal({
                cursorBlink: true,
                theme: {
                    background: '#000000',
                    foreground: '#e0e0e0',
                    cursor: '#7BC142',
                    selection: '#2D5A87'
                },
                fontFamily: 'Consolas, "Courier New", monospace',
                fontSize: 14,
                rows: 24,
                cols: 80
            });

            // Mount terminal to container
            const terminalContainer = document.getElementById('terminal-container');
            if (terminalContainer) {
                terminal.open(terminalContainer);
                console.log('✅ Terminal mounted to DOM');
                
                // Add welcome message
                terminal.write('🧬 BioCLI Companion Interactive Terminal 🧬\r\n');
                terminal.write('┌─────────────────────────────────────────────────┐\r\n');
                terminal.write('│ 🚀 Fully Interactive PowerShell Terminal       │\r\n');
                terminal.write('│ • Run real bioinformatics commands             │\r\n');
                terminal.write('│ • Get AI explanations automatically            │\r\n');
                terminal.write('│ • Learn while you work!                        │\r\n');
                terminal.write('└─────────────────────────────────────────────────┘\r\n');
                terminal.write('\r\n� Connecting to PowerShell...\r\n');
                
                // Start the shell process
                startShellProcess();
                
            } else {
                console.error('❌ Terminal container not found in DOM');
            }
            
        } else {
            throw new Error('Terminal class not available');
        }
        
    } catch (error) {
        console.error('❌ Failed to load terminal:', error);
        
        // Fallback: Show error message
        const terminalContainer = document.getElementById('terminal-container');
        if (terminalContainer) {
            terminalContainer.innerHTML = `
                <div style="color: #dc3545; padding: 20px; text-align: center;">
                    <h3>⚠️ Terminal Loading Issue</h3>
                    <p>xterm.js failed to load properly</p>
                    <p><strong>Error:</strong> ${error.message}</p>
                    <p>Check DevTools console for details</p>
                </div>
            `;
        }
    }
}

function setupBasicUI() {
    console.log('🎛️ Setting up basic UI interactions...');
    
    // Chat input functionality
    const userInput = document.getElementById('user-input');
    const sendButton = document.getElementById('send-button');
    
    if (userInput && sendButton) {
        sendButton.addEventListener('click', () => {
            const message = userInput.value.trim();
            if (message) {
                console.log('💬 User message:', message);
                addChatMessage('user', message);
                addChatMessage('ai', 'Hello! I can help explain bioinformatics commands. For now, this is a demo response.');
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
    const clearChatBtn = document.getElementById('clear-chat');
    if (clearChatBtn) {
        clearChatBtn.addEventListener('click', () => {
            const chatContainer = document.getElementById('chat-container');
            if (chatContainer) {
                // Keep welcome message, remove others
                const messages = chatContainer.querySelectorAll('.message');
                messages.forEach(msg => msg.remove());
            }
        });
    }
    
    // Toggle sidebar button
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

// Shell process management
function startShellProcess() {
    console.log('� Starting PowerShell process...');
    
    try {
        const { spawn } = require('child_process');
        
        // Spawn PowerShell process 
        shellProcess = spawn('powershell.exe', ['-NoLogo'], {
            cwd: process.cwd(),
            env: process.env
        });
        
        console.log('✅ PowerShell process started');
        terminal.write('✅ Connected to PowerShell!\r\n\r\n');
        
        // Handle shell output (stdout)
        shellProcess.stdout.on('data', (data) => {
            const output = data.toString();
            terminal.write(output);
        });
        
        // Handle shell errors (stderr)  
        shellProcess.stderr.on('data', (data) => {
            const error = data.toString();
            terminal.write(`\r\n❌ ${error}\r\n`);
        });
        
        // Handle shell process exit
        shellProcess.on('exit', (code) => {
            console.log(`🚪 PowerShell process exited with code: ${code}`);
            terminal.write(`\r\n\r\n🚪 Shell exited (code: ${code})\r\n`);
            terminal.write('🔄 Restart the app to reconnect.\r\n');
        });
        
        // Handle terminal input (user typing)
        terminal.onData((data) => {
            // Send user input to shell process
            if (shellProcess && !shellProcess.killed) {
                shellProcess.stdin.write(data);
            }
        });
        
        // Handle terminal key events for special commands
        terminal.onKey(({ key, domEvent }) => {
            // Ctrl+C handling
            if (domEvent.ctrlKey && domEvent.key === 'c') {
                if (shellProcess && !shellProcess.killed) {
                    shellProcess.kill('SIGINT');
                }
            }
        });
        
    } catch (error) {
        console.error('❌ Failed to start shell process:', error);
        terminal.write('❌ Failed to connect to PowerShell\r\n');
        terminal.write(`Error: ${error.message}\r\n`);
        terminal.write('📝 You can still use the AI chat on the right →\r\n');
    }
}

// Clean up on window close
window.addEventListener('beforeunload', () => {
    if (shellProcess && !shellProcess.killed) {
        console.log('🧹 Cleaning up shell process...');
        shellProcess.kill();
    }
});

console.log('�📝 Interactive renderer script loaded');