import { useState } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { Download, Eye } from "lucide-react";
import { useQuery } from "@tanstack/react-query";
import { apiService, Job } from "@/lib/api";
import { MolecularViewer } from "./MolecularViewer";
import { useToast } from "@/hooks/use-toast";

interface ResultsTabProps {
  selectedJobId?: string;
}

export function ResultsTab({ selectedJobId }: ResultsTabProps) {
  const [selectedPose, setSelectedPose] = useState<number>(1);
  const { toast } = useToast();

  const { data: jobData, isLoading } = useQuery({
    queryKey: ['job-status', selectedJobId],
    queryFn: () => selectedJobId ? apiService.getJobStatus(selectedJobId) : null,
    enabled: !!selectedJobId,
  });

  const getBestScore = () => {
    if (jobData?.results?.scores?.modes && jobData.results.scores.modes.length > 0) {
      return jobData.results.scores.modes[0].affinity;
    }
    if (jobData?.results?.scores?.best_affinity) {
      return jobData.results.scores.best_affinity;
    }
    return null;
  };

  const handleDownload = (format: 'pdbqt' | 'pdb') => {
    if (selectedJobId) {
      const url = apiService.getResultUrl(selectedJobId, format);
      window.open(url, '_blank');
    }
  };

  const handleViewPose = (poseNumber: number) => {
    setSelectedPose(poseNumber);
  };

  return (
    <div className="grid lg:grid-cols-3 gap-6">
      <div className="lg:col-span-2">
        <Card className="shadow-card border-border/50">
          <CardHeader>
            <CardTitle className="text-2xl bg-gradient-primary bg-clip-text text-transparent">3D Visualization</CardTitle>
            <CardDescription>Interactive molecular structure viewer</CardDescription>
          </CardHeader>
          <CardContent>
            <MolecularViewer
              proteinUrl={selectedJobId ? apiService.getDownloadUrl(selectedJobId, 'protein_prepared.pdbqt') : undefined}
              ligandUrl={selectedJobId ? apiService.getDownloadUrl(selectedJobId, 'top_pose.pdb') : undefined}
            />
          </CardContent>
        </Card>
      </div>
      
      <div className="space-y-4">
        <Card className="shadow-card border-border/50 bg-gradient-card">
          <CardHeader>
            <CardTitle className="text-xl">Docking Results</CardTitle>
            <CardDescription>Binding scores and poses</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {!selectedJobId ? (
              <div className="text-center py-8 text-muted-foreground">
                <p>Select a completed job to view results</p>
              </div>
            ) : isLoading ? (
              <div className="text-center py-8 text-muted-foreground">
                <p>Loading results...</p>
              </div>
            ) : jobData?.status !== 'completed' ? (
              <div className="text-center py-8 text-muted-foreground">
                <p>Job is not completed yet</p>
              </div>
            ) : (
              <>
                <div className="rounded-lg bg-primary/5 p-4 border border-primary/20">
                  <p className="text-sm text-muted-foreground mb-1">Best Binding Score</p>
                  <p className="text-3xl font-bold text-primary">
                    {getBestScore()?.toFixed(2) || 'N/A'}
                  </p>
                  <p className="text-xs text-muted-foreground mt-1">kcal/mol</p>
                </div>
                
                {jobData.results?.scores?.modes && jobData.results.scores.modes.length > 0 && (
                  <div>
                    <h3 className="font-semibold mb-3 text-sm text-foreground">Binding Scores</h3>
                    <div className="border rounded-lg overflow-hidden border-border/50">
                      <Table>
                        <TableHeader>
                          <TableRow className="bg-muted/30">
                            <TableHead className="text-xs">Pose</TableHead>
                            <TableHead className="text-xs">Score</TableHead>
                            <TableHead className="text-xs">RMSD</TableHead>
                            <TableHead className="text-xs">View</TableHead>
                          </TableRow>
                        </TableHeader>
                        <TableBody>
                          {jobData.results.scores.modes.map((mode, index) => (
                            <TableRow key={mode.mode} className="hover:bg-muted/20">
                              <TableCell className="font-medium text-sm">{mode.mode}</TableCell>
                              <TableCell className="text-sm">{mode.affinity.toFixed(2)}</TableCell>
                              <TableCell className="text-sm">{mode.rmsd_lb.toFixed(2)}</TableCell>
                              <TableCell>
                                <Button 
                                  variant="ghost" 
                                  size="sm" 
                                  className="h-7 w-7 p-0"
                                  onClick={() => handleViewPose(mode.mode)}
                                >
                                  <Eye className="w-3 h-3" />
                                </Button>
                              </TableCell>
                            </TableRow>
                          ))}
                        </TableBody>
                      </Table>
                    </div>
                  </div>
                )}
                
                <div className="flex gap-2 pt-2">
                  <Button 
                    variant="outline" 
                    size="sm" 
                    className="flex-1 gap-2"
                    onClick={() => handleDownload('pdbqt')}
                  >
                    <Download className="w-4 h-4" />
                    PDBQT
                  </Button>
                  <Button 
                    variant="outline" 
                    size="sm" 
                    className="flex-1 gap-2"
                    onClick={() => handleDownload('pdb')}
                  >
                    <Download className="w-4 h-4" />
                    PDB
                  </Button>
                </div>
              </>
            )}
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
