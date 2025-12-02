import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Card } from "@/components/ui/card";

interface MoleculeBuilderProps {
    onMoleculeChange: (smiles: string) => void;
}

const MoleculeBuilder = ({ onMoleculeChange }: MoleculeBuilderProps) => {
    const [selectedElement, setSelectedElement] = useState("C");
    const [smiles, setSmiles] = useState("");

    const periodicTable = [
        { symbol: "H", name: "Hydrogen", color: "bg-gray-500" },
        { symbol: "C", name: "Carbon", color: "bg-slate-600" },
        { symbol: "N", name: "Nitrogen", color: "bg-blue-600" },
        { symbol: "O", name: "Oxygen", color: "bg-red-600" },
        { symbol: "F", name: "Fluorine", color: "bg-green-500" },
        { symbol: "P", name: "Phosphorus", color: "bg-orange-500" },
        { symbol: "S", name: "Sulfur", color: "bg-yellow-600" },
        { symbol: "Cl", name: "Chlorine", color: "bg-green-600" },
        { symbol: "Br", name: "Bromine", color: "bg-red-700" },
        { symbol: "I", name: "Iodine", color: "bg-purple-600" },
    ];

    const commonMolecules = [
        { name: "Methane", smiles: "C", formula: "CH₄" },
        { name: "Ethanol", smiles: "CCO", formula: "C₂H₅OH" },
        { name: "Benzene", smiles: "c1ccccc1", formula: "C₆H₆" },
        { name: "Acetic Acid", smiles: "CC(=O)O", formula: "CH₃COOH" },
        { name: "Glucose", smiles: "C(C1C(C(C(C(O1)O)O)O)O)O", formula: "C₆H₁₂O₆" },
        { name: "Aspirin", smiles: "CC(=O)Oc1ccccc1C(=O)O", formula: "C₉H₈O₄" },
    ];

    const handleLoadTemplate = (template: string) => {
        setSmiles(template);
        onMoleculeChange(template);
    };

    return (
        <div className="h-full p-6 overflow-y-auto">
            <div className="max-w-4xl mx-auto space-y-6">
                {/* SMILES Editor */}
                <Card className="p-6 bg-black/20 border-white/10">
                    <h3 className="text-lg font-semibold text-white mb-4">SMILES Editor</h3>
                    <div className="space-y-4">
                        <div>
                            <Label className="text-white text-sm mb-2">Enter SMILES String</Label>
                            <div className="flex gap-2">
                                <Input
                                    value={smiles}
                                    onChange={(e) => setSmiles(e.target.value)}
                                    placeholder="e.g., CCO for ethanol"
                                    className="bg-white/5 border-white/10 text-white flex-1"
                                />
                                <Button onClick={() => onMoleculeChange(smiles)}>
                                    Load Molecule
                                </Button>
                            </div>
                        </div>

                        <div className="text-sm text-gray-400">
                            <p className="font-semibold text-white mb-2">SMILES Quick Guide:</p>
                            <ul className="space-y-1 text-xs">
                                <li>• C = Carbon, N = Nitrogen, O = Oxygen</li>
                                <li>• CC = Ethane, CCC = Propane</li>
                                <li>• c1ccccc1 = Benzene ring</li>
                                <li>• C=C = Double bond, C#C = Triple bond</li>
                                <li>• CCO = Ethanol (CH₃CH₂OH)</li>
                            </ul>
                        </div>
                    </div>
                </Card>

                {/* Periodic Table */}
                <Card className="p-6 bg-black/20 border-white/10">
                    <h3 className="text-lg font-semibold text-white mb-4">Periodic Table</h3>
                    <div className="grid grid-cols-5 md:grid-cols-10 gap-2">
                        {periodicTable.map((element) => (
                            <Button
                                key={element.symbol}
                                variant={selectedElement === element.symbol ? "default" : "outline"}
                                className={`h-16 flex flex-col items-center justify-center ${selectedElement === element.symbol ? element.color : "border-white/20"
                                    }`}
                                onClick={() => setSelectedElement(element.symbol)}
                            >
                                <span className="text-lg font-bold">{element.symbol}</span>
                                <span className="text-[8px] opacity-70">{element.name}</span>
                            </Button>
                        ))}
                    </div>
                    <p className="text-xs text-gray-400 mt-4">
                        Selected: <span className="text-white font-semibold">{selectedElement}</span>
                    </p>
                </Card>

                {/* Common Molecule Templates */}
                <Card className="p-6 bg-black/20 border-white/10">
                    <h3 className="text-lg font-semibold text-white mb-4">Common Molecules</h3>
                    <div className="grid grid-cols-2 md:grid-cols-3 gap-3">
                        {commonMolecules.map((mol) => (
                            <Button
                                key={mol.smiles}
                                variant="outline"
                                className="h-auto py-3 flex flex-col items-start border-purple-500/30 hover:border-purple-500/60"
                                onClick={() => handleLoadTemplate(mol.smiles)}
                            >
                                <span className="font-semibold text-white">{mol.name}</span>
                                <span className="text-xs text-gray-400">{mol.formula}</span>
                                <span className="text-[10px] text-gray-500 font-mono mt-1">
                                    {mol.smiles.substring(0, 20)}
                                    {mol.smiles.length > 20 ? "..." : ""}
                                </span>
                            </Button>
                        ))}
                    </div>
                </Card>

                {/* Instructions */}
                <Card className="p-6 bg-gradient-to-br from-cyan-500/10 to-purple-500/10 border-cyan-500/30">
                    <h3 className="text-lg font-semibold text-cyan-400 mb-3">💡 How to Build Molecules</h3>
                    <ol className="space-y-2 text-sm text-gray-300">
                        <li>1. Select an element from the periodic table</li>
                        <li>2. Enter the SMILES notation manually</li>
                        <li>3. Or choose from common molecule templates</li>
                        <li>4. Click "Load Molecule" to view in 3D</li>
                    </ol>
                    <p className="text-xs text-gray-400 mt-4">
                        Advanced 2D drawing editor coming in Phase 2!
                    </p>
                </Card>
            </div>
        </div>
    );
};

export default MoleculeBuilder;
