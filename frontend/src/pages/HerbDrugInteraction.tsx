import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Leaf, Pill, AlertTriangle, ArrowLeft, Loader2, Check, X } from "lucide-react";
import { useNavigate } from "react-router-dom";
import { useToast } from "@/hooks/use-toast";

interface HDIResults {
    herb_info: { name: string; smiles: string };
    drug_info: { name: string; smiles: string };
    interaction_risk: string;
    affected_enzymes: Array<{
        name: string;
        herb_affinity: number;
        drug_affinity: number;
        risk_level: string;
    }>;
    clinical_significance: string;
    mechanism: string;
    recommendation: string;
}

const HerbDrugInteraction = () => {
    const navigate = useNavigate();
    const { toast } = useToast();

    // Form state
    const [herbName, setHerbName] = useState("");
    const [herbCompound, setHerbCompound] = useState("");
    const [herbSmiles, setHerbSmiles] = useState("");
    const [drugName, setDrugName] = useState("");
    const [drugCompound, setDrugCompound] = useState("");
    const [drugSmiles, setDrugSmiles] = useState("");

    // Analysis state
    const [isAnalyzing, setIsAnalyzing] = useState(false);
    const [jobId, setJobId] = useState<string | null>(null);
    const [results, setResults] = useState<HDIResults | null>(null);
    const [error, setError] = useState<string | null>(null);

    const handleAnalyze = async () => {
        // Validation
        if (!herbSmiles || !drugSmiles) {
            toast({
                title: "Missing Information",
                description: "Please provide SMILES notation for both herb and drug",
                variant: "destructive"
            });
            return;
        }

        setIsAnalyzing(true);
        setError(null);
        setResults(null);

        try {
            // Prepare form data
            const formData = new FormData();
            formData.append("herb_name", herbName || "Unknown Herb");
            formData.append("herb_compound", herbCompound);
            formData.append("herb_smiles", herbSmiles);
            formData.append("drug_name", drugName || "Unknown Drug");
            formData.append("drug_compound", drugCompound);
            formData.append("drug_smiles", drugSmiles);

            // Submit job
            const submitResponse = await fetch("/api/hdi/analyze", {
                method: "POST",
                body: formData
            });

            if (!submitResponse.ok) {
                throw new Error("Failed to submit analysis");
            }

            const submitData = await submitResponse.json();
            const newJobId = submitData.job_id;
            setJobId(newJobId);

            // Poll for results
            const pollResults = async () => {
                const statusResponse = await fetch(`/api/hdi/status/${newJobId}`);
                const statusData = await statusResponse.json();

                if (statusData.status === 'completed') {
                    // Get results
                    const resultsResponse = await fetch(`/api/hdi/results/${newJobId}`);
                    const resultsData = await resultsResponse.json();
                    setResults(resultsData);
                    setIsAnalyzing(false);

                    toast({
                        title: "Analysis Complete",
                        description: "Herb-drug interaction analysis completed successfully"
                    });
                } else if (statusData.status === 'failed') {
                    throw new Error(statusData.error || "Analysis failed");
                } else {
                    // Still running, poll again
                    setTimeout(pollResults, 2000);
                }
            };

            // Start polling
            setTimeout(pollResults, 2000);

        } catch (err: any) {
            setError(err.message);
            setIsAnalyzing(false);
            toast({
                title: "Analysis Failed",
                description: err.message,
                variant: "destructive"
            });
        }
    };

    const getRiskColor = (risk: string) => {
        switch (risk) {
            case "HIGH": return "bg-red-500/20 text-red-300 border-red-500/50";
            case "MEDIUM": return "bg-yellow-500/20 text-yellow-300 border-yellow-500/50";
            case "LOW": return "bg-green-500/20 text-green-300 border-green-500/50";
            default: return "bg-gray-500/20 text-gray-300 border-gray-500/50";
        }
    };

    return (
        <div className="min-h-screen bg-background text-foreground selection:bg-primary/30">
            {/* Subtle Background */}
            <div className="fixed inset-0 pointer-events-none z-0">
                <div className="absolute top-[10%] left-[5%] w-[500px] h-[500px] bg-gradient-to-br from-green-500/20 via-emerald-500/15 to-transparent rounded-full blur-[100px] animate-pulse" />
                <div className="absolute bottom-[10%] right-[5%] w-[400px] h-[400px] bg-gradient-to-tl from-teal-500/15 via-cyan-500/10 to-transparent rounded-full blur-[90px] animate-pulse" style={{ animationDelay: '1s' }} />
            </div>

            {/* Header */}
            <header className="relative z-10 border-b border-white/10 bg-card/30 backdrop-blur-xl">
                <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6">
                    <div className="flex items-center justify-between">
                        <div className="flex items-center gap-4">
                            <Button
                                variant="ghost"
                                size="sm"
                                onClick={() => navigate("/")}
                                className="hover:bg-white/10"
                            >
                                <ArrowLeft className="w-4 h-4 mr-2" />
                                Back to Home
                            </Button>
                            <div className="h-8 w-px bg-white/20" />
                            <div className="flex items-center gap-3">
                                <div className="p-2 bg-green-500/20 rounded-lg">
                                    <Leaf className="w-6 h-6 text-green-400" />
                                </div>
                                <div>
                                    <h1 className="text-2xl font-bold text-white">Herb-Drug Interaction Analysis</h1>
                                    <p className="text-sm text-gray-400">Predict interactions between herbal compounds and pharmaceutical drugs</p>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </header>

            {/* Main Content */}
            <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-12 relative z-10">
                <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
                    {/* Input Section */}
                    <div className="lg:col-span-2 space-y-6">
                        <Card className="border-white/10 bg-card/50 backdrop-blur-xl">
                            <CardHeader>
                                <CardTitle className="flex items-center gap-2">
                                    <Leaf className="w-5 h-5 text-green-400" />
                                    Herbal Compound Input
                                </CardTitle>
                                <CardDescription>
                                    Enter the herbal compound or natural product for analysis
                                </CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-4">
                                <div className="space-y-2">
                                    <Label htmlFor="herb-name">Herb/Plant Name</Label>
                                    <Input
                                        id="herb-name"
                                        placeholder="e.g., St. John's Wort, Ginkgo biloba, Ginseng"
                                        className="bg-black/20 border-white/10"
                                        value={herbName}
                                        onChange={(e) => setHerbName(e.target.value)}
                                    />
                                </div>
                                <div className="space-y-2">
                                    <Label htmlFor="herb-compound">Active Compound (Optional)</Label>
                                    <Input
                                        id="herb-compound"
                                        placeholder="e.g., Hypericin, Ginkgolide, Ginsenoside"
                                        className="bg-black/20 border-white/10"
                                        value={herbCompound}
                                        onChange={(e) => setHerbCompound(e.target.value)}
                                    />
                                </div>
                                <div className="space-y-2">
                                    <Label htmlFor="herb-smiles">SMILES *Required</Label>
                                    <Input
                                        id="herb-smiles"
                                        placeholder="Enter SMILES notation"
                                        className="bg-black/20 border-white/10"
                                        value={herbSmiles}
                                        onChange={(e) => setHerbSmiles(e.target.value)}
                                    />
                                </div>
                            </CardContent>
                        </Card>

                        <Card className="border-white/10 bg-card/50 backdrop-blur-xl">
                            <CardHeader>
                                <CardTitle className="flex items-center gap-2">
                                    <Pill className="w-5 h-5 text-blue-400" />
                                    Pharmaceutical Drug Input
                                </CardTitle>
                                <CardDescription>
                                    Enter the pharmaceutical drug to check for interactions
                                </CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-4">
                                <div className="space-y-2">
                                    <Label htmlFor="drug-name">Drug Name</Label>
                                    <Input
                                        id="drug-name"
                                        placeholder="e.g., Warfarin, Simvastatin, Metformin"
                                        className="bg-black/20 border-white/10"
                                        value={drugName}
                                        onChange={(e) => setDrugName(e.target.value)}
                                    />
                                </div>
                                <div className="space-y-2">
                                    <Label htmlFor="drug-compound">Active Ingredient (Optional)</Label>
                                    <Input
                                        id="drug-compound"
                                        placeholder="Enter active pharmaceutical ingredient"
                                        className="bg-black/20 border-white/10"
                                        value={drugCompound}
                                        onChange={(e) => setDrugCompound(e.target.value)}
                                    />
                                </div>
                                <div className="space-y-2">
                                    <Label htmlFor="drug-smiles">SMILES *Required</Label>
                                    <Input
                                        id="drug-smiles"
                                        placeholder="Enter SMILES notation"
                                        className="bg-black/20 border-white/10"
                                        value={drugSmiles}
                                        onChange={(e) => setDrugSmiles(e.target.value)}
                                    />
                                </div>
                            </CardContent>
                        </Card>

                        <Button
                            className="w-full h-14 text-lg bg-gradient-to-r from-green-500 to-emerald-600 hover:from-green-600 hover:to-emerald-700 text-white shadow-lg"
                            onClick={handleAnalyze}
                            disabled={isAnalyzing}
                        >
                            {isAnalyzing ? (
                                <>
                                    <Loader2 className="w-5 h-5 mr-2 animate-spin" />
                                    Analyzing...
                                </>
                            ) : (
                                <>
                                    <AlertTriangle className="w-5 h-5 mr-2" />
                                    Analyze Herb-Drug Interaction
                                </>
                            )}
                        </Button>

                        {/* Results Section */}
                        {results && (
                            <Card className="border-white/10 bg-card/50 backdrop-blur-xl">
                                <CardHeader>
                                    <CardTitle className="flex items-center gap-3">
                                        <Check className="w-6 h-6 text-green-400" />
                                        Analysis Results
                                    </CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-6">
                                    {/* Risk Level */}
                                    <div>
                                        <Label className="text-sm text-gray-400 mb-2 block">Interaction Risk</Label>
                                        <Badge className={`text-lg px-4 py-2 ${getRiskColor(results.interaction_risk)}`}>
                                            {results.interaction_risk} RISK
                                        </Badge>
                                    </div>

                                    {/* Affected Enzymes */}
                                    {results.affected_enzymes && results.affected_enzymes.length > 0 && (
                                        <div>
                                            <Label className="text-sm text-gray-400 mb-3 block">Affected Enzymes</Label>
                                            <div className="space-y-2">
                                                {results.affected_enzymes.map((enzyme, idx) => (
                                                    <div key={idx} className="p-3 bg-white/5 rounded-lg">
                                                        <div className="flex items-center justify-between">
                                                            <span className="font-semibold">{enzyme.name}</span>
                                                            <Badge className={`text-xs ${getRiskColor(enzyme.risk_level)}`}>
                                                                {enzyme.risk_level}
                                                            </Badge>
                                                        </div>
                                                        <div className="text-sm text-gray-400 mt-1">
                                                            Herb: {enzyme.herb_affinity.toFixed(2)} kcal/mol |
                                                            Drug: {enzyme.drug_affinity.toFixed(2)} kcal/mol
                                                        </div>
                                                    </div>
                                                ))}
                                            </div>
                                        </div>
                                    )}

                                    {/* Clinical Significance */}
                                    <div>
                                        <Label className="text-sm text-gray-400 mb-2 block">Clinical Significance</Label>
                                        <p className="text-sm text-gray-300 whitespace-pre-wrap">{results.clinical_significance}</p>
                                    </div>

                                    {/* Mechanism */}
                                    <div>
                                        <Label className="text-sm text-gray-400 mb-2 block">Mechanism</Label>
                                        <p className="text-sm text-gray-300">{results.mechanism}</p>
                                    </div>

                                    {/* Recommendation */}
                                    <div className="bg-amber-500/10 border border-amber-500/30 rounded-lg p-4">
                                        <Label className="text-sm text-amber-400 mb-2 block flex items-center gap-2">
                                            <AlertTriangle className="w-4 h-4" />
                                            Recommendation
                                        </Label>
                                        <p className="text-sm text-gray-300">{results.recommendation}</p>
                                    </div>
                                </CardContent>
                            </Card>
                        )}

                        {/* Error Display */}
                        {error && (
                            <Card className="border-red-500/30 bg-red-500/10 backdrop-blur-xl">
                                <CardContent className="pt-6">
                                    <div className="flex items-start gap-3">
                                        <X className="w-5 h-5 text-red-400 mt-0.5" />
                                        <div>
                                            <p className="font-semibold text-red-400">Analysis Failed</p>
                                            <p className="text-sm text-gray-300 mt-1">{error}</p>
                                        </div>
                                    </div>
                                </CardContent>
                            </Card>
                        )}
                    </div>

                    {/* Info Sidebar */}
                    <div className="space-y-6">
                        <Card className="border-white/10 bg-card/50 backdrop-blur-xl">
                            <CardHeader>
                                <CardTitle className="text-lg">Analysis Targets</CardTitle>
                            </CardHeader>
                            <CardContent className="space-y-3">
                                <div className="space-y-2">
                                    <h4 className="font-semibold text-sm text-green-400">Metabolic Enzymes</h4>
                                    <ul className="text-sm text-gray-300 space-y-1">
                                        <li>• CYP3A4 (50% of drugs)</li>
                                        <li>• CYP2D6 (25% of drugs)</li>
                                        <li>• CYP2C9, CYP2C19, CYP1A2</li>
                                    </ul>
                                </div>
                                <div className="space-y-2">
                                    <h4 className="font-semibold text-sm text-blue-400">Drug Transporters</h4>
                                    <ul className="text-sm text-gray-300 space-y-1">
                                        <li>• P-glycoprotein (MDR1)</li>
                                        <li>• OATP transporters</li>
                                        <li>• BCRP, MRP proteins</li>
                                    </ul>
                                </div>
                            </CardContent>
                        </Card>

                        <Card className="border-amber-500/20 bg-amber-500/5 backdrop-blur-xl">
                            <CardHeader>
                                <CardTitle className="text-lg flex items-center gap-2 text-amber-400">
                                    <AlertTriangle className="w-5 h-5" />
                                    Common Interactions
                                </CardTitle>
                            </CardHeader>
                            <CardContent className="space-y-2 text-sm text-gray-300">
                                <p><strong className="text-amber-300">St. John's Wort:</strong> Reduces effectiveness of birth control, antidepressants</p>
                                <p><strong className="text-amber-300">Ginkgo biloba:</strong> Increases bleeding risk with anticoagulants</p>
                                <p><strong className="text-amber-300">Ginseng:</strong> May interfere with diabetes medications</p>
                                <p><strong className="text-amber-300">Garlic:</strong> Can enhance blood-thinning effects</p>
                            </CardContent>
                        </Card>

                        <Card className="border-white/10 bg-card/50 backdrop-blur-xl">
                            <CardHeader>
                                <CardTitle className="text-lg">How It Works</CardTitle>
                            </CardHeader>
                            <CardContent className="text-sm text-gray-300 space-y-2">
                                <p>1. Molecular docking against key enzymes (CYP450)</p>
                                <p>2. Binding competition analysis with drug targets</p>
                                <p>3. Transporter interference prediction</p>
                                <p>4. Risk assessment based on binding affinity</p>
                            </CardContent>
                        </Card>
                    </div>
                </div>
            </main>
        </div>
    );
};

export default HerbDrugInteraction;
