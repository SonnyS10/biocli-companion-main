const { app, BrowserWindow, Menu } = require('electron');
const path = require('path');

let mainWindow;

function createWindow() {
  console.log('🚀 Creating main window...');
  
  // Create the browser window
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 800,
    minHeight: 600,
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false
    },
    // icon: path.join(__dirname, 'assets/icon.png'), // TODO: Add custom icon
    titleBarStyle: 'default',
    show: false // Don't show until ready
  });
  
  console.log('✅ Window created, loading HTML...');

  // Load the minimal HTML file first to test basic setup
  mainWindow.loadFile('renderer/minimal.html').then(() => {
    console.log('📄 HTML file loaded successfully');
  }).catch((err) => {
    console.error('❌ Failed to load HTML file:', err);
  });

  // Show window immediately for debugging
  mainWindow.show();
  console.log('✅ Window shown');

  mainWindow.on('closed', () => {
    console.log('🚪 Window closed');
    mainWindow = null;
  });

  // Handle window closed
  mainWindow.on('closed', () => {
    mainWindow = null;
  });

  // Create application menu
  createMenu();
}

function createMenu() {
  const template = [
    {
      label: 'File',
      submenu: [
        {
          label: 'New Terminal',
          accelerator: 'Ctrl+Shift+T',
          click: () => {
            // TODO: Create new terminal tab
            mainWindow.webContents.send('new-terminal');
          }
        },
        { type: 'separator' },
        {
          label: 'Exit',
          accelerator: process.platform === 'darwin' ? 'Cmd+Q' : 'Ctrl+Q',
          click: () => {
            app.quit();
          }
        }
      ]
    },
    {
      label: 'View',
      submenu: [
        {
          label: 'Toggle AI Sidebar',
          accelerator: 'Ctrl+Shift+A',
          click: () => {
            mainWindow.webContents.send('toggle-sidebar');
          }
        },
        { type: 'separator' },
        {
          label: 'Zoom In',
          accelerator: 'Ctrl+Plus',
          click: () => {
            mainWindow.webContents.send('zoom-in');
          }
        },
        {
          label: 'Zoom Out',
          accelerator: 'Ctrl+-',
          click: () => {
            mainWindow.webContents.send('zoom-out');
          }
        },
        {
          label: 'Reset Zoom',
          accelerator: 'Ctrl+0',
          click: () => {
            mainWindow.webContents.send('zoom-reset');
          }
        }
      ]
    },
    {
      label: 'Help',
      submenu: [
        {
          label: 'About BioCLI Companion',
          click: () => {
            // TODO: Show about dialog
            mainWindow.webContents.send('show-about');
          }
        },
        {
          label: 'Documentation',
          click: () => {
            require('electron').shell.openExternal('https://github.com/SonnyS10/biocli-companion');
          }
        }
      ]
    }
  ];

  const menu = Menu.buildFromTemplate(template);
  Menu.setApplicationMenu(menu);
}

// App event listeners
app.whenReady(() => {
  console.log('🚀 Electron app ready, creating window...');
  createWindow();
});

app.on('window-all-closed', () => {
  console.log('� All windows closed');
  // On macOS, keep app running even when all windows closed
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  console.log('🔄 App activated');
  // On macOS, re-create window when dock icon clicked
  if (BrowserWindow.getAllWindows().length === 0) {
    createWindow();
  }
});