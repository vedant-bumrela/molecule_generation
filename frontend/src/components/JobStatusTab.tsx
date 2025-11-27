import { useState } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { RefreshCw, Clock, CheckCircle2, AlertCircle, Loader2, Eye } from "lucide-react";
import { useQuery } from "@tanstack/react-query";
import { apiService, Job } from "@/lib/api";
import { useToast } from "@/hooks/use-toast";


const statusConfig = {
  completed: { icon: CheckCircle2, color: "bg-success/10 text-success border-success/20", label: "Completed" },
  running: { icon: Loader2, color: "bg-secondary/10 text-secondary border-secondary/20", label: "Running" },
  pending: { icon: Clock, color: "bg-muted/10 text-muted-foreground border-muted/20", label: "Pending" },
  failed: { icon: AlertCircle, color: "bg-destructive/10 text-destructive border-destructive/20", label: "Failed" },
};

interface JobStatusTabProps {
  onViewResults?: (jobId: string) => void;
}

export function JobStatusTab({ onViewResults }: JobStatusTabProps) {
  const { toast } = useToast();

  const { data: jobsData, isLoading, refetch, isRefetching, error } = useQuery({
    queryKey: ['jobs'],
    queryFn: () => apiService.getJobs(),
    refetchInterval: 5000, // Auto-refresh every 5 seconds
  });

  // Handle errors
  if (error) {
    toast({
      title: "Error Loading Jobs",
      description: error.message,
      variant: "destructive",
    });
  }

  const handleRefresh = () => {
    refetch();
  };

  const handleViewResults = (jobId: string) => {
    if (onViewResults) {
      onViewResults(jobId);
    }
  };

  const formatTimestamp = (dateString?: string) => {
    if (!dateString) return 'Unknown';
    const date = new Date(dateString);
    const now = new Date();
    const diffMs = now.getTime() - date.getTime();
    const diffMins = Math.floor(diffMs / 60000);
    
    if (diffMins < 1) return 'Just now';
    if (diffMins < 60) return `${diffMins} minutes ago`;
    const diffHours = Math.floor(diffMins / 60);
    if (diffHours < 24) return `${diffHours} hours ago`;
    const diffDays = Math.floor(diffHours / 24);
    return `${diffDays} days ago`;
  };

  const getBestScore = (job: Job) => {
    if (job.results?.scores && job.results.scores.length > 0) {
      const bestScore = Math.min(...job.results.scores);
      return `${bestScore.toFixed(2)} kcal/mol`;
    }
    return null;
  };

  return (
    <Card className="shadow-card border-border/50">
      <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-4">
        <CardTitle className="text-2xl bg-gradient-primary bg-clip-text text-transparent">Job Status</CardTitle>
        <Button
          variant="outline"
          size="sm"
          onClick={handleRefresh}
          disabled={isRefetching}
          className="gap-2"
        >
          <RefreshCw className={`w-4 h-4 ${isRefetching ? "animate-spin" : ""}`} />
          Refresh
        </Button>
      </CardHeader>
      <CardContent>
        {isLoading ? (
          <div className="flex items-center justify-center py-8">
            <Loader2 className="w-6 h-6 animate-spin mr-2" />
            <span>Loading jobs...</span>
          </div>
        ) : jobsData?.jobs?.length === 0 ? (
          <div className="text-center py-8 text-muted-foreground">
            <p>No jobs found. Submit a job to get started.</p>
          </div>
        ) : (
          <div className="grid gap-3">
            {jobsData?.jobs?.map((job) => {
              const status = statusConfig[job.status as keyof typeof statusConfig];
              const StatusIcon = status.icon;
              const bestScore = getBestScore(job);
            
              return (
                <Card
                  key={job.job_id}
                  className="cursor-pointer hover:shadow-soft transition-all duration-300 hover:-translate-y-1 border-border/50 bg-gradient-card"
                  onClick={() => job.status === 'completed' && handleViewResults(job.job_id)}
                >
                  <CardContent className="p-4">
                    <div className="flex items-center justify-between gap-4">
                      <div className="flex items-center gap-3 flex-1">
                        <div className={`p-2 rounded-lg ${status.color.split(' ')[0]}`}>
                          <StatusIcon className={`w-5 h-5 ${job.status === "running" ? "animate-spin" : ""}`} />
                        </div>
                        <div className="flex-1 min-w-0">
                          <h3 className="font-semibold text-foreground truncate">{job.job_id}</h3>
                          <p className="text-sm text-muted-foreground">{job.job_type.toUpperCase()}</p>
                        </div>
                      </div>
                      
                      <div className="flex items-center gap-3">
                        {bestScore && (
                          <div className="text-right hidden sm:block">
                            <p className="text-sm font-semibold text-foreground">{bestScore}</p>
                            <p className="text-xs text-muted-foreground">Best Score</p>
                          </div>
                        )}
                        
                        <div className="flex flex-col items-end gap-2">
                          <div className="flex items-center gap-2">
                            <Badge className={`${status.color} border`}>
                              {status.label}
                            </Badge>
                            {job.status === 'completed' && (
                              <Button
                                variant="outline"
                                size="sm"
                                onClick={(e) => {
                                  e.stopPropagation();
                                  handleViewResults(job.job_id);
                                }}
                                className="h-6 px-2"
                              >
                                <Eye className="w-3 h-3 mr-1" />
                                View
                              </Button>
                            )}
                          </div>
                          <span className="text-xs text-muted-foreground">
                            {formatTimestamp(job.created_at)}
                          </span>
                        </div>
                      </div>
                    </div>
                    {job.message && (
                      <div className="mt-2 text-xs text-muted-foreground">
                        {job.message}
                      </div>
                    )}
                  </CardContent>
                </Card>
              );
            })}
          </div>
        )}
      </CardContent>
    </Card>
  );
}
