// BioCLI Companion Desktop App - Renderer Process
const { ipcRenderer } = require('electron');

// Terminal and UI State
let terminal;
let fitAddon;
let chatContainer;
let userInput;
let sendButton;
let sidebarVisible = true;
let currentCommand = '';

// Initialize the app when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('🚀 BioCLI Companion Desktop App Starting...');
    
    initializeUI();
    initializeTerminal();
    setupEventHandlers();
    setupIPC();
    
    console.log('✅ App initialized successfully!');
});

// Initialize UI elements
function initializeUI() {
    chatContainer = document.getElementById('chat-container');
    userInput = document.getElementById('user-input');
    sendButton = document.getElementById('send-button');
    
    // Update status
    updateStatus('ai-status', '🤖 AI Ready');
    updateStatus('terminal-status', '⚡ Terminal Ready');
    updateStatus('connection-status', '🟢 Connected');
}

// Initialize xterm.js terminal
function initializeTerminal() {
    try {
        // Create terminal instance with options
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
            lineHeight: 1.2,
            allowTransparency: true
        });

        // Create fit addon for responsive sizing
        fitAddon = new FitAddon();
        terminal.loadAddon(fitAddon);

        // Mount terminal to container
        const terminalContainer = document.getElementById('terminal-container');
        terminal.open(terminalContainer);
        
        // Fit terminal to container
        fitAddon.fit();

        // Welcome message
        showWelcomeMessage();

        // Monitor terminal input for command detection
        setupCommandMonitoring();

        console.log('✅ Terminal initialized successfully');
        updateStatus('terminal-status', '✅ Terminal Active');
        
    } catch (error) {
        console.error('❌ Failed to initialize terminal:', error);
        updateStatus('terminal-status', '❌ Terminal Error');
        
        // Fallback: Show error message
        const terminalContainer = document.getElementById('terminal-container');
        terminalContainer.innerHTML = `
            <div style="color: #dc3545; padding: 20px; text-align: center;">
                <h3>Terminal Error</h3>
                <p>Failed to initialize terminal: ${error.message}</p>
                <p>Please check console for details.</p>
            </div>
        `;
    }
}

// Show welcome message in terminal
function showWelcomeMessage() {
    const welcomeText = [
        '\r\n🧬 Welcome to BioCLI Companion Desktop! 🧬\r\n',
        '┌─────────────────────────────────────────────────┐\r\n',
        '│  Revolutionary Terminal + AI Assistant          │\r\n', 
        '│  • Run bioinformatics commands normally         │\r\n',
        '│  • Get instant AI explanations in sidebar       │\r\n',
        '│  • Learn while you work!                        │\r\n',
        '└─────────────────────────────────────────────────┘\r\n',
        '\r\n💡 Try running: bwa mem, samtools view, or fastqc\r\n',
        '🤖 AI explanations will appear automatically!\r\n\r\n'
    ].join('');

    terminal.write(welcomeText);
}

// Setup command monitoring for AI integration
function setupCommandMonitoring() {
    let currentLine = '';
    
    terminal.onKey(({ key, domEvent }) => {
        const code = domEvent.keyCode;
        
        if (code === 13) { // Enter key
            if (currentLine.trim()) {
                handleCommandEntered(currentLine.trim());
                currentLine = '';
            }
        } else if (code === 8) { // Backspace
            if (currentLine.length > 0) {
                currentLine = currentLine.slice(0, -1);
            }
        } else if (key.length === 1) { // Regular character
            currentLine += key;
        }
    });
}

// Handle when a command is entered in terminal
async function handleCommandEntered(command) {
    console.log('🔍 Command detected:', command);
    currentCommand = command;
    
    // Check if it's a bioinformatics command
    if (isBioinformaticsCommand(command)) {
        console.log('🧬 Bioinformatics command detected, requesting AI explanation...');
        await requestAIExplanation(command);
    }
}

// Check if command is bioinformatics-related
function isBioinformaticsCommand(command) {
    const bioTools = [
        'bwa', 'samtools', 'fastqc', 'star', 'bcftools', 'bedtools',
        'bowtie2', 'tophat', 'cufflinks', 'gatk', 'picard', 'trimmomatic',
        'blast', 'blastn', 'blastp', 'blastx', 'tblastn', 'tblastx',
        'muscle', 'clustalw', 'mafft', 'hisat2', 'kallisto', 'salmon'
    ];
    
    const firstWord = command.split(' ')[0].toLowerCase();
    return bioTools.includes(firstWord);
}

// Request AI explanation from backend
async function requestAIExplanation(command) {
    try {
        updateStatus('ai-status', '🤖 Thinking...');
        
        // Add user message to chat
        addChatMessage('user', `Explain: ${command}`);
        
        // Show loading in chat
        const loadingId = addLoadingMessage();
        
        // Make API request to backend (adjust URL as needed)
        const response = await axios.post('http://localhost:8000/api/explain', {
            command: command
        });
        
        // Remove loading message
        removeLoadingMessage(loadingId);
        
        // Add AI response to chat
        addChatMessage('ai', response.data.explanation);
        
        updateStatus('ai-status', '🤖 AI Ready');
        
    } catch (error) {
        console.error('❌ Failed to get AI explanation:', error);
        updateStatus('ai-status', '❌ AI Error');
        
        // Remove loading if it exists
        const loadingElements = document.querySelectorAll('.loading-message');
        loadingElements.forEach(el => el.remove());
        
        // Add error message
        addChatMessage('ai', 'Sorry, I couldn\'t explain that command right now. Please check if the backend server is running.');
    }
}

// Add message to chat interface
function addChatMessage(sender, content) {
    const messageDiv = document.createElement('div');
    messageDiv.className = `message ${sender}`;
    
    const time = new Date().toLocaleTimeString();
    
    messageDiv.innerHTML = `
        <div class="message-header">
            <span class="message-sender ${sender}">${sender === 'user' ? 'You' : 'AI Assistant'}</span>
            <span class="message-time">${time}</span>
        </div>
        <div class="message-content">${formatMessageContent(content)}</div>
    `;
    
    chatContainer.appendChild(messageDiv);
    chatContainer.scrollTop = chatContainer.scrollHeight;
}

// Format message content (handle code blocks, etc.)
function formatMessageContent(content) {
    // Basic markdown-like formatting
    return content
        .replace(/`([^`]+)`/g, '<code>$1</code>')
        .replace(/\n/g, '<br>')
        .replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>');
}

// Add loading message
function addLoadingMessage() {
    const loadingDiv = document.createElement('div');
    const loadingId = `loading-${Date.now()}`;
    loadingDiv.id = loadingId;
    loadingDiv.className = 'message ai loading-message';
    
    loadingDiv.innerHTML = `
        <div class="message-header">
            <span class="message-sender ai">AI Assistant</span>
            <span class="message-time">thinking...</span>
        </div>
        <div class="message-content">
            <span class="loading"></span> Analyzing command...
        </div>
    `;
    
    chatContainer.appendChild(loadingDiv);
    chatContainer.scrollTop = chatContainer.scrollHeight;
    
    return loadingId;
}

// Remove loading message
function removeLoadingMessage(loadingId) {
    const loadingElement = document.getElementById(loadingId);
    if (loadingElement) {
        loadingElement.remove();
    }
}

// Setup event handlers
function setupEventHandlers() {
    // Send button click
    sendButton.addEventListener('click', handleSendMessage);
    
    // Enter key in input
    userInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
            handleSendMessage();
        }
    });
    
    // Clear chat button
    document.getElementById('clear-chat').addEventListener('click', () => {
        clearChat();
    });
    
    // Clear terminal button  
    document.getElementById('clear-terminal').addEventListener('click', () => {
        if (terminal) {
            terminal.clear();
            showWelcomeMessage();
        }
    });
    
    // Toggle sidebar button
    document.getElementById('toggle-sidebar').addEventListener('click', () => {
        toggleSidebar();
    });
    
    // Window resize handler
    window.addEventListener('resize', () => {
        if (terminal && fitAddon) {
            fitAddon.fit();
        }
    });
    
    // Setup resizer functionality
    setupResizer();
}

// Handle manual message sending
async function handleSendMessage() {
    const message = userInput.value.trim();
    if (!message) return;
    
    // Clear input
    userInput.value = '';
    
    // Request AI explanation
    await requestAIExplanation(message);
}

// Clear chat messages
function clearChat() {
    const messages = chatContainer.querySelectorAll('.message');
    messages.forEach(msg => msg.remove());
    
    // Keep welcome message
    const welcomeMessage = document.querySelector('.welcome-message');
    if (!welcomeMessage) {
        addWelcomeMessage();
    }
}

// Add welcome message to chat
function addWelcomeMessage() {
    const welcomeDiv = document.createElement('div');
    welcomeDiv.className = 'welcome-message';
    welcomeDiv.innerHTML = `
        <div class="welcome-icon">🧬</div>
        <h3>Welcome to BioCLI Companion!</h3>
        <p>I'm here to help you understand bioinformatics commands. Try running a command in the terminal, and I'll automatically explain it!</p>
        <div class="example-commands">
            <p><strong>Try these examples:</strong></p>
            <code>bwa mem -t 4 reference.fa reads.fq</code><br>
            <code>samtools view -b -S alignment.sam</code><br>
            <code>fastqc *.fastq</code>
        </div>
    `;
    chatContainer.appendChild(welcomeDiv);
}

// Toggle sidebar visibility
function toggleSidebar() {
    const aiPanel = document.getElementById('ai-panel');
    const resizer = document.getElementById('resizer');
    
    sidebarVisible = !sidebarVisible;
    
    if (sidebarVisible) {
        aiPanel.classList.remove('hidden');
        resizer.classList.remove('hidden');
    } else {
        aiPanel.classList.add('hidden');
        resizer.classList.add('hidden');
    }
    
    // Resize terminal
    setTimeout(() => {
        if (terminal && fitAddon) {
            fitAddon.fit();
        }
    }, 100);
}

// Setup panel resizer
function setupResizer() {
    const resizer = document.getElementById('resizer');
    const terminalPanel = document.getElementById('terminal-panel');
    const aiPanel = document.getElementById('ai-panel');
    
    let isResizing = false;
    
    resizer.addEventListener('mousedown', (e) => {
        isResizing = true;
        resizer.classList.add('dragging');
        document.addEventListener('mousemove', handleResize);
        document.addEventListener('mouseup', stopResize);
        e.preventDefault();
    });
    
    function handleResize(e) {
        if (!isResizing) return;
        
        const containerRect = document.getElementById('app-container').getBoundingClientRect();
        const newTerminalWidth = e.clientX - containerRect.left;
        const minTerminalWidth = 400;
        const minAIWidth = 300;
        const maxTerminalWidth = containerRect.width - minAIWidth - 4; // Account for resizer width
        
        if (newTerminalWidth >= minTerminalWidth && newTerminalWidth <= maxTerminalWidth) {
            terminalPanel.style.flex = 'none';
            terminalPanel.style.width = `${newTerminalWidth}px`;
            
            // Resize terminal
            setTimeout(() => {
                if (terminal && fitAddon) {
                    fitAddon.fit();
                }
            }, 10);
        }
    }
    
    function stopResize() {
        isResizing = false;
        resizer.classList.remove('dragging');
        document.removeEventListener('mousemove', handleResize);
        document.removeEventListener('mouseup', stopResize);
    }
}

// Setup IPC communication with main process
function setupIPC() {
    // Handle menu commands
    ipcRenderer.on('toggle-sidebar', () => {
        toggleSidebar();
    });
    
    ipcRenderer.on('new-terminal', () => {
        // TODO: Implement multiple terminal tabs
        console.log('New terminal requested');
    });
    
    ipcRenderer.on('zoom-in', () => {
        changeZoom(0.1);
    });
    
    ipcRenderer.on('zoom-out', () => {
        changeZoom(-0.1);
    });
    
    ipcRenderer.on('zoom-reset', () => {
        setZoom(1.0);
    });
}

// Zoom functionality
let currentZoom = 1.0;

function changeZoom(delta) {
    currentZoom = Math.max(0.5, Math.min(2.0, currentZoom + delta));
    setZoom(currentZoom);
}

function setZoom(zoom) {
    currentZoom = zoom;
    document.body.style.zoom = zoom;
    updateStatus('zoom-level', `${Math.round(zoom * 100)}%`);
    
    // Resize terminal after zoom
    setTimeout(() => {
        if (terminal && fitAddon) {
            fitAddon.fit();
        }
    }, 100);
}

// Update status bar
function updateStatus(elementId, text) {
    const element = document.getElementById(elementId);
    if (element) {
        element.textContent = text;
    }
}