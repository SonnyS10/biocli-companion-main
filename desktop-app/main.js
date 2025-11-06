const { app, BrowserWindow } = require('electron');
const path = require('path');

console.log('🚀 BioCLI Companion starting...');

let mainWindow;

function createWindow() {
  console.log('📱 Creating window...');
  
  mainWindow = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false
    }
  });

  console.log('📄 Loading our beautiful app...');
  mainWindow.loadFile('renderer/index.html');
  
  console.log('🪟 Window created and loading...');
  
  mainWindow.on('closed', () => {
    console.log('🚪 Window closed');
    mainWindow = null;
  });
}

app.whenReady().then(() => {
  console.log('⚡ App ready!');
  createWindow();
});

app.on('window-all-closed', () => {
  console.log('🚪 All windows closed');
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

console.log('📝 Main script loaded');