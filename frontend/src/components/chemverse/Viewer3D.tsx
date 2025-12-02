import { useEffect, useRef, useState } from "react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Card } from "@/components/ui/card";
import { Maximize2, RotateCw, ZoomIn, ZoomOut, Download } from "lucide-react";

interface Viewer3DProps {
    molecule: string | null;
}

const Viewer3D = ({ molecule }: Viewer3DProps) => {
    const viewerRef = useRef<HTMLDivElement>(null);
    const [viewerInstance, setViewerInstance] = useState<any>(null);
    const [smiles, setSmiles] = useState(molecule || "");
    const [style, setStyle] = useState<"stick" | "sphere" | "cartoon">("stick");
    const [loading, setLoading] = useState(false);

    // Load 3Dmol.js library dynamically
    useEffect(() => {
        const script = document.createElement("script");
        script.src = "https://3Dmol.csb.pitt.edu/build/3Dmol-min.js";
        script.async = true;
        document.body.appendChild(script);

        return () => {
            document.body.removeChild(script);
        };
    }, []);

    // Initialize viewer
    useEffect(() => {
        if (!viewerRef.current || !window.$3Dmol) return;

        const config = {
            backgroundColor: "0x0a0e27",
        };

        const viewer = window.$3Dmol.createViewer(viewerRef.current, config);
        setViewerInstance(viewer);

        return () => {
            viewer.clear();
        };
    }, [viewerRef.current, window.$3Dmol]);

    // Load molecule when SMILES changes
    useEffect(() => {
        if (molecule && molecule !== smiles) {
            setSmiles(molecule);
            loadMolecule(molecule);
        }
    }, [molecule]);

    const loadMolecule = async (smilesInput: string) => {
        if (!viewerInstance || !smilesInput) return;

        setLoading(true);
        try {
            // Convert SMILES to 3D structure via backend
            const response = await fetch("/api/chemverse/smiles-to-3d", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({ smiles: smilesInput }),
            });

            if (!response.ok) throw new Error("Failed to convert SMILES");

            const data = await response.json();
            const pdbData = data.pdb;

            // Clear previous molecule
            viewerInstance.clear();

            // Add new molecule
            viewerInstance.addModel(pdbData, "pdb");

            // Apply style
            const styleConfig = getStyleConfig(style);
            viewerInstance.setStyle({}, styleConfig);

            // Add labels for atoms (optional)
            viewerInstance.addPropertyLabels("elem", {}, { fontColor: "white", fontSize: 12 });

            // Render
            viewerInstance.zoomTo();
            viewerInstance.render();
        } catch (error) {
            console.error("Error loading molecule:", error);
        } finally {
            setLoading(false);
        }
    };

    const getStyleConfig = (styleType: string) => {
        switch (styleType) {
            case "stick":
                return { stick: { colorscheme: "Jmol" } };
            case "sphere":
                return { sphere: { colorscheme: "Jmol" } };
            case "cartoon":
                return { cartoon: { color: "spectrum" } };
            default:
                return { stick: {} };
        }
    };

    const handleStyleChange = (newStyle: "stick" | "sphere" | "cartoon") => {
        setStyle(newStyle);
        if (viewerInstance) {
            const styleConfig = getStyleConfig(newStyle);
            viewerInstance.setStyle({}, styleConfig);
            viewerInstance.render();
        }
    };

    const handleZoomIn = () => {
        if (viewerInstance) {
            viewerInstance.zoom(1.2);
            viewerInstance.render();
        }
    };

    const handleZoomOut = () => {
        if (viewerInstance) {
            viewerInstance.zoom(0.8);
            viewerInstance.render();
        }
    };

    const handleReset = () => {
        if (viewerInstance) {
            viewerInstance.zoomTo();
            viewerInstance.render();
        }
    };

    const handleScreenshot = () => {
        if (viewerInstance) {
            const imgData = viewerInstance.pngURI();
            const link = document.createElement("a");
            link.href = imgData;
            link.download = "molecule.png";
            link.click();
        }
    };

    return (
        <div className="h-full flex flex-col">
            {/* Controls */}
            <div className="p-4 border-b border-white/10">
                <div className="flex gap-4 items-end">
                    <div className="flex-1">
                        <Label className="text-white text-sm mb-2">SMILES Input</Label>
                        <Input
                            value={smiles}
                            onChange={(e) => setSmiles(e.target.value)}
                            placeholder="Enter SMILES (e.g., c1ccccc1)"
                            className="bg-white/5 border-white/10 text-white"
                            onKeyPress={(e) => {
                                if (e.key === "Enter") loadMolecule(smiles);
                            }}
                        />
                    </div>
                    <Button onClick={() => loadMolecule(smiles)} disabled={loading}>
                        {loading ? "Loading..." : "Load"}
                    </Button>
                </div>

                <div className="flex gap-2 mt-4">
                    <div className="flex gap-1">
                        <Button
                            size="sm"
                            variant={style === "stick" ? "default" : "outline"}
                            onClick={() => handleStyleChange("stick")}
                            className="text-xs"
                        >
                            Stick
                        </Button>
                        <Button
                            size="sm"
                            variant={style === "sphere" ? "default" : "outline"}
                            onClick={() => handleStyleChange("sphere")}
                            className="text-xs"
                        >
                            Sphere
                        </Button>
                    </div>

                    <div className="flex gap-1 ml-auto">
                        <Button size="icon" variant="outline" onClick={handleZoomIn}>
                            <ZoomIn className="w-4 h-4" />
                        </Button>
                        <Button size="icon" variant="outline" onClick={handleZoomOut}>
                            <ZoomOut className="w-4 h-4" />
                        </Button>
                        <Button size="icon" variant="outline" onClick={handleReset}>
                            <RotateCw className="w-4 h-4" />
                        </Button>
                        <Button size="icon" variant="outline" onClick={handleScreenshot}>
                            <Download className="w-4 h-4" />
                        </Button>
                    </div>
                </div>
            </div>

            {/* 3D Viewer Canvas */}
            <div className="flex-1 relative">
                <div ref={viewerRef} className="w-full h-full" />
                {!smiles && (
                    <div className="absolute inset-0 flex items-center justify-center">
                        <Card className="p-6 bg-black/60 backdrop-blur border-white/10">
                            <p className="text-gray-400 text-center">
                                Enter a SMILES string or select a quick example
                                <br />
                                <span className="text-xs text-gray-500 mt-2 block">
                                    Try: c1ccccc1 (benzene) or CC(=O)O (acetic acid)
                                </span>
                            </p>
                        </Card>
                    </div>
                )}
            </div>
        </div>
    );
};

// Declare global $3Dmol
declare global {
    interface Window {
        $3Dmol: any;
    }
}

export default Viewer3D;
