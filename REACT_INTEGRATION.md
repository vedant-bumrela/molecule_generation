# React Frontend Integration

This document explains how the React frontend has been integrated with the Flask backend for the molecular docking system.

## Architecture

- **Frontend**: React app with TypeScript, Vite, TailwindCSS, and shadcn/ui components
- **Backend**: Flask API serving both the React build files and API endpoints
- **3D Visualization**: NGL.js for molecular structure visualization
- **State Management**: TanStack Query for API data fetching and caching

## Quick Start

### Option 1: Integrated Startup (Recommended)
```powershell
# Build frontend and start integrated server
.\start_integrated.ps1
```

### Option 2: Manual Steps
```powershell
# 1. Build the React frontend
.\build_frontend.ps1

# 2. Start the Flask server
python -m src.api.app
```

## Access Points

- **Web Application**: http://localhost:5000
- **API Endpoints**: http://localhost:5000/api/*

## Features

### Frontend Components

1. **SubmitJobTab**: Submit docking jobs with file uploads and advanced options
2. **JobStatusTab**: View and monitor job progress with real-time updates
3. **ResultsTab**: Interactive 3D visualization and results analysis
4. **BatchProcessingTab**: Handle multiple compounds at once
5. **MolecularViewer**: 3D molecular structure visualization using NGL.js

### API Integration

- Real-time job status updates every 5 seconds
- File upload handling for protein and ligand structures
- SMILES string support for ligand input
- Download functionality for result files
- Error handling and user feedback

### Key Technologies

- **React 18** with TypeScript
- **Vite** for fast development and building
- **TailwindCSS** for styling
- **shadcn/ui** for UI components
- **TanStack Query** for API state management
- **NGL.js** for 3D molecular visualization
- **Lucide React** for icons

## Development

### Frontend Development
```powershell
cd frontend
npm install
npm run dev  # Development server on http://localhost:5173
```

### Backend Development
```powershell
python -m src.api.app  # Flask server on http://localhost:5000
```

## File Structure

```
frontend/
├── src/
│   ├── components/          # React components
│   │   ├── ui/             # shadcn/ui components
│   │   ├── SubmitJobTab.tsx
│   │   ├── JobStatusTab.tsx
│   │   ├── ResultsTab.tsx
│   │   ├── BatchProcessingTab.tsx
│   │   └── MolecularViewer.tsx
│   ├── lib/
│   │   ├── api.ts          # API service layer
│   │   └── utils.ts
│   ├── pages/
│   │   └── Index.tsx       # Main application page
│   └── main.tsx           # Application entry point
├── dist/                  # Built files (served by Flask)
└── package.json
```

## Environment Variables

Create a `.env` file in the `frontend/` directory:

```
VITE_API_BASE_URL=http://localhost:5000/api
```

## Troubleshooting

### Frontend Build Issues
- Ensure Node.js is installed (v18+ recommended)
- Delete `node_modules` and run `npm install` again
- Check for TypeScript errors in the console

### API Connection Issues
- Verify Flask server is running on port 5000
- Check CORS settings in Flask app
- Ensure API endpoints are accessible at `/api/*`

### 3D Visualization Issues
- NGL.js loads from CDN - ensure internet connection
- Check browser console for JavaScript errors
- Verify molecular files are in correct format (PDB/PDBQT)

## Next Steps

1. **Performance Optimization**: Implement code splitting and lazy loading
2. **Enhanced Visualization**: Add more NGL.js features and controls
3. **Real-time Updates**: WebSocket integration for live job progress
4. **Testing**: Add unit and integration tests
5. **Docker Integration**: Update Dockerfile to include React build process
