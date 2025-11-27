import { useEffect, useRef, useState } from 'react';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Eye, EyeOff, RotateCcw } from 'lucide-react';

// Import NGL
declare global {
  interface Window {
    NGL: any;
  }
}

interface MolecularViewerProps {
  proteinUrl?: string;
  ligandUrl?: string;
  className?: string;
}

export function MolecularViewer({ proteinUrl, ligandUrl, className = '' }: MolecularViewerProps) {
  const viewerRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<any>(null);
  const proteinComponentRef = useRef<any>(null);
  const ligandComponentRef = useRef<any>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [showProtein, setShowProtein] = useState(true);
  const [showLigand, setShowLigand] = useState(true);
  const [representation, setRepresentation] = useState<'cartoon' | 'surface'>('cartoon');

  useEffect(() => {
    // Load NGL dynamically
    const loadNGL = async () => {
      if (window.NGL) return;
      
      const script = document.createElement('script');
      script.src = 'https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.js';
      script.async = true;
      
      return new Promise<void>((resolve, reject) => {
        script.onload = () => resolve();
        script.onerror = reject;
        document.head.appendChild(script);
      });
    };

    const initViewer = async () => {
      try {
        await loadNGL();
        
        if (viewerRef.current && window.NGL) {
          // Initialize NGL Stage
          stageRef.current = new window.NGL.Stage(viewerRef.current, {
            backgroundColor: 'white',
            cameraType: 'perspective',
            cameraFov: 40
          });

          // Handle resize
          const handleResize = () => {
            if (stageRef.current) {
              stageRef.current.handleResize();
            }
          };

          window.addEventListener('resize', handleResize);
          
          return () => {
            window.removeEventListener('resize', handleResize);
          };
        }
      } catch (error) {
        console.error('Failed to load NGL:', error);
      }
    };

    initViewer();
  }, []);

  useEffect(() => {
    if (!stageRef.current) return;

    const loadStructures = async () => {
      setIsLoading(true);
      
      try {
        // Clear existing components
        stageRef.current.removeAllComponents();
        proteinComponentRef.current = null;
        ligandComponentRef.current = null;

        // Load protein
        if (proteinUrl) {
          try {
            const proteinBlob = await fetch(proteinUrl).then(r => r.blob());
            const proteinFile = new File([proteinBlob], 'protein.pdb', { type: 'chemical/x-pdb' });
            
            proteinComponentRef.current = await stageRef.current.loadFile(proteinFile, { ext: 'pdb' });
            updateProteinRepresentation();
          } catch (error) {
            console.warn('Error loading protein:', error);
          }
        }

        // Load ligand
        if (ligandUrl) {
          try {
            const ligandBlob = await fetch(ligandUrl).then(r => r.blob());
            const ligandFile = new File([ligandBlob], 'ligand.pdb', { type: 'chemical/x-pdb' });
            
            ligandComponentRef.current = await stageRef.current.loadFile(ligandFile, { ext: 'pdb' });
            ligandComponentRef.current.addRepresentation('ball+stick', {
              multipleBond: true,
              colorScheme: 'element',
              sele: '/0' // Select first model (best pose)
            });
          } catch (error) {
            console.warn('Error loading ligand:', error);
          }
        }

        // Auto-view
        if (proteinComponentRef.current || ligandComponentRef.current) {
          stageRef.current.autoView();
        }
      } catch (error) {
        console.error('Error loading structures:', error);
      } finally {
        setIsLoading(false);
      }
    };

    if (proteinUrl || ligandUrl) {
      loadStructures();
    }
  }, [proteinUrl, ligandUrl, representation]);

  const updateProteinRepresentation = () => {
    if (!proteinComponentRef.current) return;

    proteinComponentRef.current.removeAllRepresentations();

    if (representation === 'cartoon') {
      proteinComponentRef.current.addRepresentation('cartoon', { color: 'chainname' });
      proteinComponentRef.current.addRepresentation('licorice', {
        sele: 'hetero and not water',
        multipleBond: true
      });
    } else if (representation === 'surface') {
      proteinComponentRef.current.addRepresentation('surface', {
        opacity: 0.7,
        colorScheme: 'bfactor'
      });
    }
  };

  const toggleProteinVisibility = () => {
    if (proteinComponentRef.current) {
      proteinComponentRef.current.setVisibility(!showProtein);
      setShowProtein(!showProtein);
    }
  };

  const toggleLigandVisibility = () => {
    if (ligandComponentRef.current) {
      ligandComponentRef.current.setVisibility(!showLigand);
      setShowLigand(!showLigand);
    }
  };

  const resetView = () => {
    if (stageRef.current) {
      stageRef.current.autoView();
    }
  };

  const changeRepresentation = (newRep: 'cartoon' | 'surface') => {
    setRepresentation(newRep);
  };

  return (
    <Card className={`shadow-card hover:shadow-soft transition-shadow duration-300 ${className}`}>
      <CardHeader>
        <CardTitle className="text-xl">3D Molecular Viewer</CardTitle>
      </CardHeader>
      <CardContent>
        <div className="space-y-4">
          {/* Viewer Container */}
          <div className="relative">
            <div
              ref={viewerRef}
              className="w-full h-96 border border-border rounded-lg bg-white"
              style={{ minHeight: '400px' }}
            />
            {isLoading && (
              <div className="absolute inset-0 flex items-center justify-center bg-background/80 rounded-lg">
                <div className="flex items-center gap-2">
                  <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-primary"></div>
                  <span className="text-sm text-muted-foreground">Loading structures...</span>
                </div>
              </div>
            )}
          </div>

          {/* Controls */}
          <div className="flex flex-wrap gap-2">
            <div className="flex gap-1">
              <Button
                variant="outline"
                size="sm"
                onClick={toggleProteinVisibility}
                className="flex items-center gap-1"
              >
                {showProtein ? <Eye className="w-4 h-4" /> : <EyeOff className="w-4 h-4" />}
                Protein
              </Button>
              <Button
                variant="outline"
                size="sm"
                onClick={toggleLigandVisibility}
                className="flex items-center gap-1"
              >
                {showLigand ? <Eye className="w-4 h-4" /> : <EyeOff className="w-4 h-4" />}
                Ligand
              </Button>
            </div>

            <div className="flex gap-1">
              <Button
                variant={representation === 'cartoon' ? 'default' : 'outline'}
                size="sm"
                onClick={() => changeRepresentation('cartoon')}
              >
                Cartoon
              </Button>
              <Button
                variant={representation === 'surface' ? 'default' : 'outline'}
                size="sm"
                onClick={() => changeRepresentation('surface')}
              >
                Surface
              </Button>
            </div>

            <Button
              variant="outline"
              size="sm"
              onClick={resetView}
              className="flex items-center gap-1"
            >
              <RotateCcw className="w-4 h-4" />
              Reset View
            </Button>
          </div>

          {!proteinUrl && !ligandUrl && (
            <div className="text-center py-8 text-muted-foreground">
              <p>No structures loaded. Submit a job to view results.</p>
            </div>
          )}
        </div>
      </CardContent>
    </Card>
  );
}
