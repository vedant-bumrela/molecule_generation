/**
 * API service for molecular docking system
 * Handles all communication with Flask backend
 */

export interface Job {
  job_id: string;
  job_type: string;
  status: 'pending' | 'running' | 'completed' | 'failed';
  message?: string;
  results?: {
    scores?: {
      modes?: Array<{
        mode: number;
        affinity: number;
        rmsd_lb: number;
        rmsd_ub: number;
      }>;
      best_affinity?: number;
    };
    output_path?: string;
    top_pose_path?: string;
    num_poses?: number;
    [key: string]: any;
  };
  created_at?: string;
  updated_at?: string;
}

export interface JobSubmissionData {
  job_type: string;
  protein: File;
  ligand?: File;
  smiles?: string;
  box_center_x?: number;
  box_center_y?: number;
  box_center_z?: number;
  box_size_x?: number;
  box_size_y?: number;
  box_size_z?: number;
  exhaustiveness?: number;
  num_modes?: number;
  run_rescoring?: boolean;
  num_samples?: number;
  confidence_model?: boolean;
}

export interface BatchJobData {
  protein: File;
  library: File;
  method: string;
  num_workers: number;
  center_x?: number;
  center_y?: number;
  center_z?: number;
  size_x?: number;
  size_y?: number;
  size_z?: number;
  exhaustiveness?: number;
  num_modes?: number;
  run_rescoring?: boolean;
  num_samples?: number;
}

class ApiService {
  private baseUrl: string;

  constructor() {
    this.baseUrl = (import.meta as any).env?.VITE_API_BASE_URL || 'http://localhost:5000/api';
  }

  /**
   * Submit a docking job
   */
  async submitJob(data: JobSubmissionData): Promise<{ job_id: string }> {
    const formData = new FormData();
    
    // Add required fields
    formData.append('job_type', data.job_type);
    formData.append('protein', data.protein);
    
    // Add ligand input (either file or SMILES)
    if (data.ligand) {
      formData.append('ligand', data.ligand);
    }
    if (data.smiles) {
      formData.append('smiles', data.smiles);
    }
    
    // Add optional parameters
    if (data.box_center_x !== undefined) formData.append('box_center_x', data.box_center_x.toString());
    if (data.box_center_y !== undefined) formData.append('box_center_y', data.box_center_y.toString());
    if (data.box_center_z !== undefined) formData.append('box_center_z', data.box_center_z.toString());
    if (data.box_size_x !== undefined) formData.append('box_size_x', data.box_size_x.toString());
    if (data.box_size_y !== undefined) formData.append('box_size_y', data.box_size_y.toString());
    if (data.box_size_z !== undefined) formData.append('box_size_z', data.box_size_z.toString());
    if (data.exhaustiveness !== undefined) formData.append('exhaustiveness', data.exhaustiveness.toString());
    if (data.num_modes !== undefined) formData.append('num_modes', data.num_modes.toString());
    if (data.run_rescoring !== undefined) formData.append('run_rescoring', data.run_rescoring.toString());
    if (data.num_samples !== undefined) formData.append('num_samples', data.num_samples.toString());
    if (data.confidence_model !== undefined) formData.append('confidence_model', data.confidence_model.toString());

    const response = await fetch(`${this.baseUrl}/submit`, {
      method: 'POST',
      body: formData,
    });

    if (!response.ok) {
      const error = await response.json().catch(() => ({ error: 'Unknown error' }));
      throw new Error(error.error || `HTTP error ${response.status}`);
    }

    return response.json();
  }

  /**
   * Submit a batch job
   */
  async submitBatchJob(data: BatchJobData): Promise<{ job_id: string }> {
    const formData = new FormData();
    
    formData.append('protein', data.protein);
    formData.append('library', data.library);
    formData.append('method', data.method);
    formData.append('num_workers', data.num_workers.toString());
    
    // Add optional parameters
    if (data.center_x !== undefined) formData.append('center_x', data.center_x.toString());
    if (data.center_y !== undefined) formData.append('center_y', data.center_y.toString());
    if (data.center_z !== undefined) formData.append('center_z', data.center_z.toString());
    if (data.size_x !== undefined) formData.append('size_x', data.size_x.toString());
    if (data.size_y !== undefined) formData.append('size_y', data.size_y.toString());
    if (data.size_z !== undefined) formData.append('size_z', data.size_z.toString());
    if (data.exhaustiveness !== undefined) formData.append('exhaustiveness', data.exhaustiveness.toString());
    if (data.num_modes !== undefined) formData.append('num_modes', data.num_modes.toString());
    if (data.run_rescoring !== undefined) formData.append('run_rescoring', data.run_rescoring.toString());
    if (data.num_samples !== undefined) formData.append('num_samples', data.num_samples.toString());

    const response = await fetch(`${this.baseUrl}/submit/batch`, {
      method: 'POST',
      body: formData,
    });

    if (!response.ok) {
      const error = await response.json().catch(() => ({ error: 'Unknown error' }));
      throw new Error(error.error || `HTTP error ${response.status}`);
    }

    return response.json();
  }

  /**
   * Get all jobs
   */
  async getJobs(): Promise<{ jobs: Job[] }> {
    const response = await fetch(`${this.baseUrl}/jobs`);
    
    if (!response.ok) {
      throw new Error(`HTTP error ${response.status}`);
    }

    return response.json();
  }

  /**
   * Get job status
   */
  async getJobStatus(jobId: string): Promise<Job> {
    const response = await fetch(`${this.baseUrl}/status/${jobId}`);
    
    if (!response.ok) {
      throw new Error(`HTTP error ${response.status}`);
    }

    return response.json();
  }

  /**
   * Get batch job results
   */
  async getBatchResults(jobId: string): Promise<{ results: any[] }> {
    const response = await fetch(`${this.baseUrl}/batch/results/${jobId}`);
    
    if (!response.ok) {
      throw new Error(`HTTP error ${response.status}`);
    }

    return response.json();
  }

  /**
   * Download file from job results
   */
  getDownloadUrl(jobId: string, filename: string): string {
    return `${this.baseUrl}/download/${jobId}/${filename}`;
  }

  /**
   * Get result file URL
   */
  getResultUrl(jobId: string, format: 'pdbqt' | 'pdb'): string {
    return `${this.baseUrl}/result/${jobId}/${format}`;
  }

  /**
   * Download batch results CSV
   */
  getBatchResultsDownloadUrl(jobId: string): string {
    return `${this.baseUrl}/download/${jobId}/docking_results.csv`;
  }
}

export const apiService = new ApiService();
