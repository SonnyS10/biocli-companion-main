const { app, BrowserWindow } = require('electron');

console.log('🚀 Starting minimal Electron app...');

app.whenReady().then(() => {
  console.log('📱 App ready, creating window...');
  
  const win = new BrowserWindow({
    width: 800,
    height: 600,
    show: true // Show immediately
  });
  
  console.log('🪟 Window created');
  
  win.loadURL('data:text/html;charset=utf-8,<h1>SUCCESS!</h1><p>Minimal Electron app is working!</p>');
  
  console.log('✅ App fully loaded');
});

app.on('window-all-closed', () => {
  console.log('🚪 All windows closed');
  // Don't quit on macOS
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  console.log('🔄 App activated');
});

console.log('📝 Main.js script loaded');