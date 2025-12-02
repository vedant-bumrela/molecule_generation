import { useState } from "react";
import { useNavigate } from "react-router-dom";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Beaker, Atom, FlaskConical, Lightbulb, Home } from "lucide-react";
import Viewer3D from "@/components/chemverse/Viewer3D";
import MoleculeBuilder from "@/components/chemverse/MoleculeBuilder";
import ReactionLab from "@/components/chemverse/ReactionLab";

const ChemVerse = () => {
    const navigate = useNavigate();
    const [activeTab, setActiveTab] = useState("viewer");
    const [currentMolecule, setCurrentMolecule] = useState<string | null>(null);

    return (
        <div className="min-h-screen bg-gradient-to-br from-slate-950 via-slate-900 to-indigo-950">
            {/* Header */}
            <nav className="border-b border-white/10 bg-black/20 backdrop-blur-xl sticky top-0 z-50">
                <div className="container mx-auto px-6 py-4">
                    <div className="flex items-center justify-between">
                        <div className="flex items-center gap-3">
                            <div className="w-10 h-10 rounded-lg bg-gradient-to-br from-cyan-500 to-purple-600 flex items-center justify-center">
                                <Atom className="w-6 h-6 text-white" />
                            </div>
                            <div>
                                <h1 className="text-2xl font-bold text-white">ChemVerse</h1>
                                <p className="text-xs text-gray-400">Virtual Chemistry Lab</p>
                            </div>
                        </div>
                        <div className="flex gap-2">
                            <Button
                                variant="outline"
                                className="border-white/50 text-white hover:bg-white/10"
                                onClick={() => navigate("/")}
                            >
                                <Home className="w-4 h-4 mr-2" />
                                Back to Home
                            </Button>
                            <Button variant="outline" className="border-cyan-500/50 text-cyan-400 hover:bg-cyan-500/10">
                                <Lightbulb className="w-4 h-4 mr-2" />
                                Tutorial
                            </Button>
                        </div>
                    </div>
                </div>
            </nav>

            {/* Main Content */}
            <main className="container mx-auto px-6 py-8">
                <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
                    {/* Sidebar */}
                    <div className="lg:col-span-1 space-y-4">
                        <Card className="border-white/10 bg-black/40 backdrop-blur-xl">
                            <CardHeader>
                                <CardTitle className="text-lg text-white flex items-center gap-2">
                                    <Beaker className="w-5 h-5 text-cyan-400" />
                                    Lab Tools
                                </CardTitle>
                            </CardHeader>
                            <CardContent className="space-y-2">
                                <Button
                                    variant={activeTab === "viewer" ? "default" : "ghost"}
                                    className="w-full justify-start"
                                    onClick={() => setActiveTab("viewer")}
                                >
                                    <Atom className="w-4 h-4 mr-2" />
                                    3D Viewer
                                </Button>
                                <Button
                                    variant={activeTab === "builder" ? "default" : "ghost"}
                                    className="w-full justify-start"
                                    onClick={() => setActiveTab("builder")}
                                >
                                    <FlaskConical className="w-4 h-4 mr-2" />
                                    Molecule Builder
                                </Button>
                                <Button
                                    variant={activeTab === "reactor" ? "default" : "ghost"}
                                    className="w-full justify-start"
                                    onClick={() => setActiveTab("reactor")}
                                >
                                    <Beaker className="w-4 h-4 mr-2" />
                                    Reaction Lab
                                </Button>
                            </CardContent>
                        </Card>

                        {/* Quick Examples */}
                        <Card className="border-white/10 bg-black/40 backdrop-blur-xl">
                            <CardHeader>
                                <CardTitle className="text-sm text-white">Quick Start</CardTitle>
                            </CardHeader>
                            <CardContent className="space-y-2 text-sm">
                                <Button
                                    variant="outline"
                                    size="sm"
                                    className="w-full text-xs border-purple-500/30"
                                    onClick={() => setCurrentMolecule("CC(=O)Oc1ccccc1C(=O)O")}
                                >
                                    Load Aspirin
                                </Button>
                                <Button
                                    variant="outline"
                                    size="sm"
                                    className="w-full text-xs border-purple-500/30"
                                    onClick={() => setCurrentMolecule("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")}
                                >
                                    Load Caffeine
                                </Button>
                                <Button
                                    variant="outline"
                                    size="sm"
                                    className="w-full text-xs border-purple-500/30"
                                    onClick={() => setCurrentMolecule("c1ccccc1")}
                                >
                                    Load Benzene
                                </Button>
                            </CardContent>
                        </Card>
                    </div>

                    {/* Main Canvas */}
                    <div className="lg:col-span-3">
                        <Card className="border-white/10 bg-black/40 backdrop-blur-xl h-[calc(100vh-12rem)]">
                            <CardContent className="p-0 h-full">
                                <Tabs value={activeTab} onValueChange={setActiveTab} className="h-full">
                                    <TabsContent value="viewer" className="h-full m-0">
                                        <Viewer3D molecule={currentMolecule} />
                                    </TabsContent>
                                    <TabsContent value="builder" className="h-full m-0">
                                        <MoleculeBuilder onMoleculeChange={setCurrentMolecule} />
                                    </TabsContent>
                                    <TabsContent value="reactor" className="h-full m-0">
                                        <ReactionLab />
                                    </TabsContent>
                                </Tabs>
                            </CardContent>
                        </Card>
                    </div>
                </div>
            </main>
        </div>
    );
};

export default ChemVerse;
