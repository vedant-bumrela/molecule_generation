import { useState } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Input } from "@/components/ui/input";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Button } from "@/components/ui/button";
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from "@/components/ui/accordion";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { Download, Loader2 } from "lucide-react";
import { useToast } from "@/hooks/use-toast";

const mockBatchJobs = [
  { id: "batch_001", method: "Vina", status: "completed", progress: 100 },
  { id: "batch_002", method: "DiffDock", status: "running", progress: 65 },
];

export function BatchProcessingTab() {
  const [method, setMethod] = useState("vina");
  const { toast } = useToast();

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    toast({
      title: "Batch Job Submitted",
      description: "Your batch processing job has been queued.",
    });
  };

  return (
    <div className="space-y-6">
      <Card className="shadow-card border-border/50 bg-gradient-card">
        <CardHeader>
          <CardTitle className="text-2xl bg-gradient-primary bg-clip-text text-transparent">Batch Processing</CardTitle>
          <CardDescription>Process multiple compounds against a single protein target</CardDescription>
        </CardHeader>
        <CardContent>
          <form onSubmit={handleSubmit} className="space-y-6">
            <div className="grid md:grid-cols-2 gap-4">
              <div className="space-y-2">
                <Label htmlFor="batchProtein">Protein Structure (PDB)</Label>
                <Input id="batchProtein" type="file" accept=".pdb" className="bg-background" />
                <p className="text-xs text-muted-foreground">Upload a protein structure in PDB format</p>
              </div>
              
              <div className="space-y-2">
                <Label htmlFor="batchLibrary">Compound Library</Label>
                <Input id="batchLibrary" type="file" accept=".csv,.sdf" className="bg-background" />
                <p className="text-xs text-muted-foreground">CSV with SMILES column or SDF file</p>
              </div>
            </div>

            <div className="grid md:grid-cols-2 gap-4">
              <div className="space-y-2">
                <Label htmlFor="batchMethod">Docking Method</Label>
                <Select value={method} onValueChange={setMethod}>
                  <SelectTrigger id="batchMethod" className="bg-background">
                    <SelectValue />
                  </SelectTrigger>
                  <SelectContent className="bg-popover z-50">
                    <SelectItem value="vina">AutoDock Vina</SelectItem>
                    <SelectItem value="gnina">GNINA</SelectItem>
                    <SelectItem value="equibind">EquiBind</SelectItem>
                    <SelectItem value="diffdock">DiffDock</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              
              <div className="space-y-2">
                <Label htmlFor="numWorkers">Number of Workers</Label>
                <Input id="numWorkers" type="number" min="1" max="16" defaultValue="4" className="bg-background" />
                <p className="text-xs text-muted-foreground">Parallel workers for processing</p>
              </div>
            </div>

            <Accordion type="single" collapsible className="w-full">
              <AccordionItem value="advanced" className="border-border/50">
                <AccordionTrigger className="text-sm font-semibold">Advanced Options</AccordionTrigger>
                <AccordionContent className="space-y-4 pt-4">
                  <div className="grid grid-cols-2 md:grid-cols-3 gap-3">
                    <div className="space-y-2">
                      <Label htmlFor="bCenterX" className="text-xs">Center X</Label>
                      <Input id="bCenterX" type="number" step="0.1" defaultValue="0" className="bg-background h-9" />
                    </div>
                    <div className="space-y-2">
                      <Label htmlFor="bSizeX" className="text-xs">Size X</Label>
                      <Input id="bSizeX" type="number" step="0.1" defaultValue="20" className="bg-background h-9" />
                    </div>
                    <div className="space-y-2">
                      <Label htmlFor="bExh" className="text-xs">Exhaustiveness</Label>
                      <Input id="bExh" type="number" min="1" max="32" defaultValue="8" className="bg-background h-9" />
                    </div>
                  </div>
                </AccordionContent>
              </AccordionItem>
            </Accordion>

            <Button type="submit" size="lg" className="w-full bg-gradient-primary hover:opacity-90 transition-opacity shadow-glow">
              Submit Batch Job
            </Button>
          </form>
        </CardContent>
      </Card>

      <Card className="shadow-card border-border/50">
        <CardHeader>
          <CardTitle className="text-xl">Batch Jobs</CardTitle>
          <CardDescription>Track your batch processing jobs</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="border rounded-lg overflow-hidden border-border/50">
            <Table>
              <TableHeader>
                <TableRow className="bg-muted/30">
                  <TableHead>Job ID</TableHead>
                  <TableHead>Method</TableHead>
                  <TableHead>Status</TableHead>
                  <TableHead>Progress</TableHead>
                  <TableHead>Actions</TableHead>
                </TableRow>
              </TableHeader>
              <TableBody>
                {mockBatchJobs.map((job) => (
                  <TableRow key={job.id} className="hover:bg-muted/20">
                    <TableCell className="font-medium">{job.id}</TableCell>
                    <TableCell>{job.method}</TableCell>
                    <TableCell>
                      <Badge className={job.status === "completed" ? "bg-success/10 text-success border-success/20" : "bg-secondary/10 text-secondary border-secondary/20"}>
                        {job.status === "running" && <Loader2 className="w-3 h-3 mr-1 animate-spin" />}
                        {job.status}
                      </Badge>
                    </TableCell>
                    <TableCell>
                      <div className="flex items-center gap-2">
                        <Progress value={job.progress} className="w-20" />
                        <span className="text-xs text-muted-foreground">{job.progress}%</span>
                      </div>
                    </TableCell>
                    <TableCell>
                      {job.status === "completed" && (
                        <Button variant="outline" size="sm" className="gap-2">
                          <Download className="w-3 h-3" />
                          Results
                        </Button>
                      )}
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
