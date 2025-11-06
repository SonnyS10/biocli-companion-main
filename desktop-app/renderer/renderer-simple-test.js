// BioCLI Companion - Simple Working Terminal
console.log('🚀 Loading Simple BioCLI Terminal...');

let terminal;

document.addEventListener('DOMContentLoaded', () => {
    console.log('📄 DOM loaded, creating simple functional terminal...');
    
    initializeSimpleTerminal();
    setupBasicUI();
    
    console.log('✅ Simple terminal ready!');
});

function initializeSimpleTerminal() {
    console.log('🖥️ Creating simple terminal...');
    
    try {
        if (typeof Terminal === 'undefined') {
            throw new Error('xterm.js not loaded');
        }
        
        // Create basic terminal
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
        
        // Show simple prompt
        showSimplePrompt();
        
        // Setup simple input handling
        setupSimpleInput();
        
    } catch (error) {
        console.error('❌ Terminal creation failed:', error);
        showError(error);
    }
}

function showSimplePrompt() {
    terminal.clear();
    terminal.write('🧬 BioCLI Companion - Simple Terminal 🧬\r\n');
    terminal.write('═══════════════════════════════════════\r\n');
    terminal.write('🚀 Terminal UI is working!\r\n');
    terminal.write('💡 Type "help" for available commands\r\n');
    terminal.write('🤖 AI integration coming next!\r\n\r\n');
    terminal.write('biocli> ');
}

function setupSimpleInput() {
    let currentCommand = '';
    
    terminal.onData((data) => {
        const char = data;
        
        if (char === '\r') {
            // Enter pressed
            terminal.write('\r\n');
            
            if (currentCommand.trim()) {
                handleCommand(currentCommand.trim());
            }
            
            currentCommand = '';
            terminal.write('biocli> ');
            
        } else if (char === '\x7f' || char === '\b') {
            // Backspace
            if (currentCommand.length > 0) {
                currentCommand = currentCommand.slice(0, -1);
                terminal.write('\b \b');
            }
            
        } else if (char >= ' ') {
            // Regular character
            currentCommand += char;
            terminal.write(char);
        }
    });
}

function handleCommand(command) {
    console.log('📝 Command entered:', command);
    
    // Simple built-in commands
    switch (command.toLowerCase()) {
        case 'help':
            terminal.write('Available commands:\r\n');
            terminal.write('  help     - Show this help\r\n');
            terminal.write('  clear    - Clear terminal\r\n');
            terminal.write('  echo <text> - Echo text\r\n');
            terminal.write('  bwa      - Test bioinformatics detection\r\n');
            terminal.write('  samtools - Test bioinformatics detection\r\n');
            terminal.write('  fastqc   - Test bioinformatics detection\r\n');
            break;
            
        case 'clear':
            showSimplePrompt();
            return; // Don't show prompt again
            
        case 'bwa':
        case 'samtools':
        case 'fastqc':
        case 'star':
        case 'bcftools':
            terminal.write(`🧬 Bioinformatics command detected: ${command}\r\n`);
            addBioCommandToChat(command);
            break;
            
        default:
            if (command.startsWith('echo ')) {
                const text = command.substring(5);
                terminal.write(`${text}\r\n`);
            } else {
                terminal.write(`Command not recognized: ${command}\r\n`);
                terminal.write('Type "help" for available commands\r\n');
            }
    }
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