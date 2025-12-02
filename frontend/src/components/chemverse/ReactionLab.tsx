import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Card } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Beaker, ArrowRight, Zap, AlertTriangle } from "lucide-react";
import { useToast } from "@/hooks/use-toast";

const ReactionLab = () => {
    const { toast } = useToast();
    const [reactant1, setReactant1] = useState("");
    const [reactant2, setReactant2] = useState("");
    const [reactionType, setReactionType] = useState("");
    const [products, setProducts] = useState<string[]>([]);
    const [isReacting, setIsReacting] = useState(false);
    const [reactionData, setReactionData] = useState<any>(null);

    const reactionTemplates = [
        { value: "acid_base", label: "Acid-Base Neutralization", example: "HCl + NaOH → NaCl + H2O" },
        { value: "esterification", label: "Ester Formation", example: "R-COOH + R'-OH → R-COO-R' + H2O" },
        { value: "addition", label: "Addition Reaction", example: "C=C + H2 → C-C-H" },
        { value: "substitution", label: "Substitution", example: "R-X + Nu → R-Nu + X" },
    ];

    const commonReactions = [
        {
            name: "Aspirin Synthesis",
            reactant1: "c1ccc(cc1)O",
            reactant2: "CC(=O)O",
            type: "esterification",
        },
        {
            name: "Water Formation",
            reactant1: "[H][H]",
            reactant2: "O=O",
            type: "addition",
        },
    ];

    const handleReact = async () => {
        if (!reactant1 || !reactant2 || !reactionType) {
            toast({
                title: "Missing Information",
                description: "Please provide both reactants and select a reaction type",
                variant: "destructive",
            });
            return;
        }

        setIsReacting(true);
        try {
            const response = await fetch("/api/chemverse/react", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({
                    reactants: [reactant1, reactant2],
                    reaction_type: reactionType,
                }),
            });

            if (!response.ok) throw new Error("Reaction failed");

            const data = await response.json();
            setProducts(data.products);
            setReactionData(data);

            toast({
                title: "Reaction Complete!",
                description: `Generated ${data.products.length} product(s)`,
            });
        } catch (error) {
            toast({
                title: "Reaction Failed",
                description: "Could not predict reaction products",
                variant: "destructive",
            });
        } finally {
            setIsReacting(false);
        }
    };

    const loadExample = (example: any) => {
        setReactant1(example.reactant1);
        setReactant2(example.reactant2);
        setReactionType(example.type);
    };

    return (
        <div className="h-full p-6 overflow-y-auto">
            <div className="max-w-6xl mx-auto space-y-6">
                {/* Reaction Setup */}
                <Card className="p-6 bg-black/20 border-white/10">
                    <div className="flex items-center gap-3 mb-6">
                        <Beaker className="w-6 h-6 text-cyan-400" />
                        <h3 className="text-xl font-semibold text-white">Reaction Lab</h3>
                    </div>

                    <div className="grid md:grid-cols-3 gap-6">
                        {/* Reactant 1 */}
                        <div>
                            <Label className="text-white text-sm mb-2">Reactant 1 (SMILES)</Label>
                            <Input
                                value={reactant1}
                                onChange={(e) => setReactant1(e.target.value)}
                                placeholder="e.g., CCO"
                                className="bg-white/5 border-white/10 text-white"
                            />
                        </div>

                        {/* Reaction Type */}
                        <div>
                            <Label className="text-white text-sm mb-2">Reaction Type</Label>
                            <Select value={reactionType} onValueChange={setReactionType}>
                                <SelectTrigger className="bg-white/5 border-white/10 text-white">
                                    <SelectValue placeholder="Select reaction" />
                                </SelectTrigger>
                                <SelectContent>
                                    {reactionTemplates.map((template) => (
                                        <SelectItem key={template.value} value={template.value}>
                                            {template.label}
                                        </SelectItem>
                                    ))}
                                </SelectContent>
                            </Select>
                        </div>

                        {/* Reactant 2 */}
                        <div>
                            <Label className="text-white text-sm mb-2">Reactant 2 (SMILES)</Label>
                            <Input
                                value={reactant2}
                                onChange={(e) => setReactant2(e.target.value)}
                                placeholder="e.g., CC(=O)O"
                                className="bg-white/5 border-white/10 text-white"
                            />
                        </div>
                    </div>

                    <Button
                        onClick={handleReact}
                        disabled={isReacting}
                        className="w-full mt-6 bg-gradient-to-r from-cyan-500 to-purple-600 hover:from-cyan-600 hover:to-purple-700"
                        size="lg"
                    >
                        {isReacting ? (
                            <>
                                <Zap className="w-4 h-4 mr-2 animate-pulse" />
                                Reacting...
                            </>
                        ) : (
                            <>
                                <Beaker className="w-4 h-4 mr-2" />
                                Perform Reaction
                            </>
                        )}
                    </Button>
                </Card>

                {/* Reaction Visualization */}
                {products.length > 0 && reactionData && (
                    <Card className="p-6 bg-black/20 border-white/10">
                        <h3 className="text-lg font-semibold text-white mb-4">Reaction Results</h3>

                        <div className="flex items-center justify-center gap-4 mb-6">
                            <div className="text-center">
                                <p className="text-xs text-gray-400 mb-2">Reactants</p>
                                <div className="space-y-2">
                                    <Badge variant="outline" className="font-mono text-xs">
                                        {reactant1}
                                    </Badge>
                                    <p className="text-gray-500 text-sm">+</p>
                                    <Badge variant="outline" className="font-mono text-xs">
                                        {reactant2}
                                    </Badge>
                                </div>
                            </div>

                            <ArrowRight className="w-8 h-8 text-cyan-400" />

                            <div className="text-center">
                                <p className="text-xs text-gray-400 mb-2">Products</p>
                                <div className="space-y-2">
                                    {products.map((product, idx) => (
                                        <Badge key={idx} variant="secondary" className="font-mono text-xs bg-purple-500/20">
                                            {product}
                                        </Badge>
                                    ))}
                                </div>
                            </div>
                        </div>

                        {/* Reaction Details */}
                        <div className="grid md:grid-cols-3 gap-4 mt-6">
                            <Card className="p-4 bg-white/5 border-white/10">
                                <p className="text-xs text-gray-400 mb-1">Feasibility</p>
                                <p className={`text-lg font-semibold ${reactionData.feasible ? "text-green-400" : "text-red-400"}`}>
                                    {reactionData.feasible ? "✓ Feasible" : "✗ Not Feasible"}
                                </p>
                            </Card>

                            <Card className="p-4 bg-white/5 border-white/10">
                                <p className="text-xs text-gray-400 mb-1">Energy Change (ΔG)</p>
                                <p className="text-lg font-semibold text-cyan-400">
                                    {reactionData.energy_change?.toFixed(2)} kcal/mol
                                </p>
                            </Card>

                            <Card className="p-4 bg-white/5 border-white/10">
                                <p className="text-xs text-gray-400 mb-1">Mechanism</p>
                                <p className="text-sm text-white">{reactionData.mechanism}</p>
                            </Card>
                        </div>

                        {reactionData.safety_warnings && reactionData.safety_warnings.length > 0 && (
                            <Card className="p-4 bg-amber-500/10 border-amber-500/30 mt-4">
                                <div className="flex items-start gap-3">
                                    <AlertTriangle className="w-5 h-5 text-amber-400 mt-0.5" />
                                    <div>
                                        <p className="font-semibold text-amber-400 mb-2">Safety Warnings</p>
                                        <ul className="space-y-1 text-sm text-gray-300">
                                            {reactionData.safety_warnings.map((warning: string, idx: number) => (
                                                <li key={idx}>• {warning}</li>
                                            ))}
                                        </ul>
                                    </div>
                                </div>
                            </Card>
                        )}
                    </Card>
                )}

                {/* Example Reactions */}
                <Card className="p-6 bg-black/20 border-white/10">
                    <h3 className="text-lg font-semibold text-white mb-4">Example Reactions</h3>
                    <div className="grid md:grid-cols-2 gap-3">
                        {commonReactions.map((example, idx) => (
                            <Button
                                key={idx}
                                variant="outline"
                                className="h-auto py-4 flex flex-col items-start border-purple-500/30"
                                onClick={() => loadExample(example)}
                            >
                                <span className="font-semibold text-white">{example.name}</span>
                                <span className="text-xs text-gray-400 mt-1">
                                    {reactionTemplates.find((t) => t.value === example.type)?.label}
                                </span>
                            </Button>
                        ))}
                    </div>
                </Card>
            </div>
        </div>
    );
};

export default ReactionLab;
