import { useState } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Input } from "@/components/ui/input";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { RadioGroup, RadioGroupItem } from "@/components/ui/radio-group";
import { Checkbox } from "@/components/ui/checkbox";
import { Button } from "@/components/ui/button";
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from "@/components/ui/accordion";
import { Beaker, Dna, Sparkles, Loader2 } from "lucide-react";
import { useToast } from "@/hooks/use-toast";
import { apiService, JobSubmissionData } from "@/lib/api";
import { useMutation } from "@tanstack/react-query";

const methodInfo = {
  vina: {
    title: "AutoDock Vina",
    description: "Classical molecular docking engine that uses a scoring function and search algorithm to predict binding poses and affinities.",
    icon: Beaker,
    features: ["Speed: Fast", "Accuracy: Good baseline", "Best for: Standard docking, virtual screening"]
  },
  equibind: {
    title: "EquiBind (Fast ML)",
    description: "Fast machine learning approach that directly predicts protein-ligand complexes without requiring iterative optimization.",
    icon: Sparkles,
    features: ["Speed: Very Fast", "Accuracy: ML-powered", "Best for: High-throughput screening"]
  },
  diffdock: {
    title: "DiffDock (Generative ML)",
    description: "State-of-the-art generative diffusion model for molecular docking with high accuracy predictions.",
    icon: Dna,
    features: ["Speed: Moderate", "Accuracy: State-of-the-art", "Best for: Precision docking, novel compounds"]
  }
};

export function SubmitJobTab() {
  const [jobType, setJobType] = useState("vina");
  const [ligandInputType, setLigandInputType] = useState("file");
  const [showAdvanced, setShowAdvanced] = useState(false);
  const { toast } = useToast();

  const submitJobMutation = useMutation({
    mutationFn: (data: JobSubmissionData) => apiService.submitJob(data),
    onSuccess: (result) => {
      toast({
        title: "Job Submitted Successfully",
        description: `Your docking job has been queued. Job ID: ${result.job_id}`,
      });
    },
    onError: (error: Error) => {
      toast({
        title: "Submission Failed",
        description: error.message,
        variant: "destructive",
      });
    },
  });

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    
    const form = e.target as HTMLFormElement;
    const formData = new FormData(form);
    
    // Get form values
    const proteinFile = (form.elements.namedItem('proteinFile') as HTMLInputElement)?.files?.[0];
    const ligandFile = (form.elements.namedItem('ligandFile') as HTMLInputElement)?.files?.[0];
    const smilesInput = (form.elements.namedItem('smilesInput') as HTMLInputElement)?.value;
    
    // Validation
    if (!proteinFile) {
      toast({
        title: "Validation Error",
        description: "Please select a protein file",
        variant: "destructive",
      });
      return;
    }
    
    if (ligandInputType === 'file' && !ligandFile) {
      toast({
        title: "Validation Error",
        description: "Please select a ligand file",
        variant: "destructive",
      });
      return;
    }
    
    if (ligandInputType === 'smiles' && !smilesInput?.trim()) {
      toast({
        title: "Validation Error",
        description: "Please enter a SMILES string",
        variant: "destructive",
      });
      return;
    }

    // Prepare submission data
    const submissionData: JobSubmissionData = {
      job_type: jobType,
      protein: proteinFile,
    };

    // Add ligand data
    if (ligandInputType === 'file' && ligandFile) {
      submissionData.ligand = ligandFile;
    } else if (ligandInputType === 'smiles' && smilesInput) {
      submissionData.smiles = smilesInput.trim();
    }

    // Add advanced options if enabled
    if (showAdvanced) {
      const getNumberValue = (name: string) => {
        const value = (form.elements.namedItem(name) as HTMLInputElement)?.value;
        return value ? parseFloat(value) : undefined;
      };
      
      const getBooleanValue = (name: string) => {
        return (form.elements.namedItem(name) as HTMLInputElement)?.checked;
      };

      submissionData.box_center_x = getNumberValue('centerX');
      submissionData.box_center_y = getNumberValue('centerY');
      submissionData.box_center_z = getNumberValue('centerZ');
      submissionData.box_size_x = getNumberValue('sizeX');
      submissionData.box_size_y = getNumberValue('sizeY');
      submissionData.box_size_z = getNumberValue('sizeZ');

      if (jobType === 'vina') {
        submissionData.exhaustiveness = getNumberValue('exhaustiveness');
        submissionData.num_modes = getNumberValue('numModes');
        submissionData.run_rescoring = getBooleanValue('rescoring');
      } else if (jobType === 'diffdock') {
        submissionData.num_samples = getNumberValue('numSamples');
        submissionData.confidence_model = getBooleanValue('confidence');
      }
    }

    submitJobMutation.mutate(submissionData);
  };

  const selectedMethod = methodInfo[jobType as keyof typeof methodInfo];
  const MethodIcon = selectedMethod.icon;

  return (
    <div className="grid lg:grid-cols-2 gap-6">
      <Card className="shadow-card hover:shadow-soft transition-shadow duration-300 bg-gradient-card border-border/50">
        <CardHeader>
          <CardTitle className="text-2xl bg-gradient-primary bg-clip-text text-transparent">Submit Docking Job</CardTitle>
          <CardDescription>Configure your molecular docking parameters</CardDescription>
        </CardHeader>
        <CardContent>
          <form onSubmit={handleSubmit} className="space-y-6">
            <div className="space-y-2">
              <Label htmlFor="jobType">Job Type</Label>
              <Select value={jobType} onValueChange={setJobType}>
                <SelectTrigger id="jobType" className="bg-background">
                  <SelectValue />
                </SelectTrigger>
                <SelectContent className="bg-popover z-50">
                  <SelectItem value="vina">AutoDock Vina (Classical)</SelectItem>
                  <SelectItem value="equibind">EquiBind (Fast ML)</SelectItem>
                  <SelectItem value="diffdock">DiffDock (Generative ML)</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <div className="space-y-2">
              <Label htmlFor="proteinFile">Protein Structure (PDB)</Label>
              <Input id="proteinFile" name="proteinFile" type="file" accept=".pdb" className="bg-background" />
            </div>

            <div className="space-y-3">
              <Label>Ligand Input</Label>
              <RadioGroup value={ligandInputType} onValueChange={setLigandInputType}>
                <div className="flex items-center space-x-2">
                  <RadioGroupItem value="file" id="ligandFile" />
                  <Label htmlFor="ligandFile" className="font-normal cursor-pointer">Upload Ligand File</Label>
                </div>
                <div className="flex items-center space-x-2">
                  <RadioGroupItem value="smiles" id="ligandSmiles" />
                  <Label htmlFor="ligandSmiles" className="font-normal cursor-pointer">Enter SMILES String</Label>
                </div>
              </RadioGroup>
            </div>

            {ligandInputType === "file" ? (
              <div className="space-y-2">
                <Label htmlFor="ligandFile">Ligand File (SDF, MOL2, PDB)</Label>
                <Input id="ligandFile" name="ligandFile" type="file" accept=".sdf,.mol2,.pdb,.pdbqt" className="bg-background" />
              </div>
            ) : (
              <div className="space-y-2">
                <Label htmlFor="smilesInput">SMILES String</Label>
                <Input id="smilesInput" name="smilesInput" placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O" className="bg-background" />
              </div>
            )}

            <div className="flex items-center space-x-2">
              <Checkbox id="advanced" checked={showAdvanced} onCheckedChange={(checked) => setShowAdvanced(checked as boolean)} />
              <Label htmlFor="advanced" className="font-normal cursor-pointer">Show Advanced Options</Label>
            </div>

            {showAdvanced && (
              <Accordion type="single" collapsible className="w-full">
                <AccordionItem value="advanced" className="border-border/50">
                  <AccordionTrigger className="text-sm font-semibold">Advanced Options</AccordionTrigger>
                  <AccordionContent className="space-y-4 pt-4">
                    <div className="grid grid-cols-2 gap-3">
                      <div className="space-y-2">
                        <Label htmlFor="centerX" className="text-xs">Center X</Label>
                        <Input id="centerX" name="centerX" type="number" step="0.1" className="bg-background h-9" />
                      </div>
                      <div className="space-y-2">
                        <Label htmlFor="sizeX" className="text-xs">Size X</Label>
                        <Input id="sizeX" name="sizeX" type="number" step="0.1" defaultValue="20" className="bg-background h-9" />
                      </div>
                      <div className="space-y-2">
                        <Label htmlFor="centerY" className="text-xs">Center Y</Label>
                        <Input id="centerY" name="centerY" type="number" step="0.1" className="bg-background h-9" />
                      </div>
                      <div className="space-y-2">
                        <Label htmlFor="sizeY" className="text-xs">Size Y</Label>
                        <Input id="sizeY" name="sizeY" type="number" step="0.1" defaultValue="20" className="bg-background h-9" />
                      </div>
                      <div className="space-y-2">
                        <Label htmlFor="centerZ" className="text-xs">Center Z</Label>
                        <Input id="centerZ" name="centerZ" type="number" step="0.1" className="bg-background h-9" />
                      </div>
                      <div className="space-y-2">
                        <Label htmlFor="sizeZ" className="text-xs">Size Z</Label>
                        <Input id="sizeZ" name="sizeZ" type="number" step="0.1" defaultValue="20" className="bg-background h-9" />
                      </div>
                    </div>

                    {jobType === "vina" && (
                      <div className="space-y-3 pt-2">
                        <div className="space-y-2">
                          <Label htmlFor="exhaustiveness" className="text-sm">Exhaustiveness</Label>
                          <Input id="exhaustiveness" name="exhaustiveness" type="number" min="1" max="32" defaultValue="8" className="bg-background" />
                        </div>
                        <div className="space-y-2">
                          <Label htmlFor="numModes" className="text-sm">Number of Poses</Label>
                          <Input id="numModes" name="numModes" type="number" min="1" max="20" defaultValue="9" className="bg-background" />
                        </div>
                        <div className="flex items-center space-x-2">
                          <Checkbox id="rescoring" name="rescoring" />
                          <Label htmlFor="rescoring" className="text-sm font-normal cursor-pointer">Run GNINA ML Rescoring</Label>
                        </div>
                      </div>
                    )}

                    {jobType === "diffdock" && (
                      <div className="space-y-3 pt-2">
                        <div className="space-y-2">
                          <Label htmlFor="numSamples" className="text-sm">Number of Samples</Label>
                          <Input id="numSamples" name="numSamples" type="number" min="1" max="20" defaultValue="10" className="bg-background" />
                        </div>
                        <div className="flex items-center space-x-2">
                          <Checkbox id="confidence" name="confidence" defaultChecked />
                          <Label htmlFor="confidence" className="text-sm font-normal cursor-pointer">Use Confidence Model</Label>
                        </div>
                      </div>
                    )}
                  </AccordionContent>
                </AccordionItem>
              </Accordion>
            )}

            <Button 
              type="submit" 
              size="lg" 
              className="w-full bg-gradient-primary hover:opacity-90 transition-opacity shadow-glow"
              disabled={submitJobMutation.isPending}
            >
              {submitJobMutation.isPending ? (
                <>
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                  Submitting...
                </>
              ) : (
                'Submit Job'
              )}
            </Button>
          </form>
        </CardContent>
      </Card>

      <Card className="shadow-card hover:shadow-soft transition-shadow duration-300 bg-gradient-card border-border/50">
        <CardHeader>
          <div className="flex items-center gap-3">
            <div className="p-3 bg-primary/10 rounded-lg">
              <MethodIcon className="w-6 h-6 text-primary" />
            </div>
            <div>
              <CardTitle className="text-xl">{selectedMethod.title}</CardTitle>
              <CardDescription>Method Information</CardDescription>
            </div>
          </div>
        </CardHeader>
        <CardContent className="space-y-4">
          <p className="text-sm text-muted-foreground leading-relaxed">{selectedMethod.description}</p>
          <ul className="space-y-2">
            {selectedMethod.features.map((feature, index) => (
              <li key={index} className="flex items-start gap-2 text-sm">
                <span className="text-primary mt-1">•</span>
                <span className="text-foreground">{feature}</span>
              </li>
            ))}
          </ul>
        </CardContent>
      </Card>
    </div>
  );
}
