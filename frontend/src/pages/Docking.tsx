import { useState } from "react";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { SubmitJobTab } from "@/components/SubmitJobTab";
import { JobStatusTab } from "@/components/JobStatusTab";
import { ResultsTab } from "@/components/ResultsTab";
import { BatchProcessingTab } from "@/components/BatchProcessingTab";
import { Atom, Activity, FileText, Layers, ArrowLeft, FlaskConical } from "lucide-react";
import { Button } from "@/components/ui/button";
import { useNavigate } from "react-router-dom";

const Docking = () => {
    const navigate = useNavigate();
    const [selectedJobId, setSelectedJobId] = useState<string>();
    const [activeTab, setActiveTab] = useState("submit");

    const handleViewResults = (jobId: string) => {
        setSelectedJobId(jobId);
        setActiveTab("results");
    };

    return (
        <div className="min-h-screen bg-background text-foreground selection:bg-primary/30">
            {/* Subtle Background */}
            <div className="fixed inset-0 pointer-events-none z-0">
                <div className="absolute top-[10%] left-[5%] w-[500px] h-[500px] bg-gradient-to-br from-primary/20 via-purple-500/15 to-transparent rounded-full blur-[100px] animate-pulse" />
                <div className="absolute bottom-[10%] right-[5%] w-[400px] h-[400px] bg-gradient-to-tl from-secondary/15 via-cyan-500/10 to-transparent rounded-full blur-[90px] animate-pulse" style={{ animationDelay: '1s' }} />
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
                                <div className="p-2 bg-primary/20 rounded-lg">
                                    <FlaskConical className="w-6 h-6 text-primary" />
                                </div>
                                <div>
                                    <h1 className="text-2xl font-bold text-white">Molecular Docking</h1>
                                    <p className="text-sm text-gray-400">AI-powered drug discovery platform</p>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </header>

            {/* Main Content */}
            <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-12 relative z-10">
                <div className="bg-card/50 backdrop-blur-xl border border-white/10 rounded-3xl shadow-2xl p-6 md:p-8 animate-fade-in">
                    <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
                        <TabsList className="grid w-full grid-cols-2 lg:grid-cols-4 h-auto gap-4 bg-black/20 p-2 rounded-2xl border border-white/5">
                            <TabsTrigger
                                value="submit"
                                className="data-[state=active]:bg-primary data-[state=active]:text-white data-[state=active]:shadow-glow transition-all duration-300 py-3 rounded-xl text-sm font-medium tracking-wide"
                            >
                                <Atom className="w-4 h-4 mr-2" />
                                Submit Job
                            </TabsTrigger>
                            <TabsTrigger
                                value="status"
                                className="data-[state=active]:bg-primary data-[state=active]:text-white data-[state=active]:shadow-glow transition-all duration-300 py-3 rounded-xl text-sm font-medium tracking-wide"
                            >
                                <Activity className="w-4 h-4 mr-2" />
                                Job Status
                            </TabsTrigger>
                            <TabsTrigger
                                value="results"
                                className="data-[state=active]:bg-primary data-[state=active]:text-white data-[state=active]:shadow-glow transition-all duration-300 py-3 rounded-xl text-sm font-medium tracking-wide"
                            >
                                <FileText className="w-4 h-4 mr-2" />
                                Results
                            </TabsTrigger>
                            <TabsTrigger
                                value="batch"
                                className="data-[state=active]:bg-primary data-[state=active]:text-white data-[state=active]:shadow-glow transition-all duration-300 py-3 rounded-xl text-sm font-medium tracking-wide"
                            >
                                <Layers className="w-4 h-4 mr-2" />
                                Batch Processing
                            </TabsTrigger>
                        </TabsList>

                        <div className="mt-8 min-h-[400px]">
                            <TabsContent value="submit" className="animate-fade-in focus-visible:outline-none focus-visible:ring-0">
                                <SubmitJobTab />
                            </TabsContent>

                            <TabsContent value="status" className="animate-fade-in focus-visible:outline-none focus-visible:ring-0">
                                <JobStatusTab onViewResults={handleViewResults} />
                            </TabsContent>

                            <TabsContent value="results" className="animate-fade-in focus-visible:outline-none focus-visible:ring-0">
                                <ResultsTab selectedJobId={selectedJobId} />
                            </TabsContent>

                            <TabsContent value="batch" className="animate-fade-in focus-visible:outline-none focus-visible:ring-0">
                                <BatchProcessingTab />
                            </TabsContent>
                        </div>
                    </Tabs>
                </div>
            </main>

            {/* Footer Decoration */}
            <div className="fixed bottom-0 left-0 right-0 h-1 bg-gradient-to-r from-transparent via-primary/50 to-transparent opacity-50" />
        </div>
    );
};

export default Docking;
