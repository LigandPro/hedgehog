// Screen types
export type Screen = 
  | 'welcome'
  | 'configMain'
  | 'configDescriptors'
  | 'configFilters'
  | 'configSynthesis'
  | 'configDocking'
  | 'pipelineRunner';

// Config types
export interface MainConfig {
  generated_mols_path: string;
  target_mols_path: string;
  folder_to_save: string;
  n_jobs: number;
  sample_size: number;
  batch_size: number;
  save_sampled_mols: boolean;
  pains_file_path: string;
  mcf_file_path: string;
  ligand_preparation_tool?: string;
  protein_preparation_tool?: string;
  config_descriptors: string;
  config_structFilters: string;
  config_synthesis: string;
  config_docking: string;
}

export interface DescriptorBorder {
  min: number | string;
  max: number | string;
}

export interface DescriptorsConfig {
  run: boolean;
  batch_size: number;
  filter_data: boolean;
  borders: Record<string, DescriptorBorder | string[] | boolean>;
  filtered_cols_to_plot: string[];
  discrete_features_to_plot: string[];
  not_to_smooth_plot_by_sides: string[];
  renamer: Record<string, string>;
}

export interface FiltersConfig {
  run: boolean;
  run_before_descriptors: boolean;
  filter_data: boolean;
  alerts_data_path: string;
  calculate_common_alerts: boolean;
  include_rulesets: string[];
  exclude_descriptions: Record<string, string[]>;
  calculate_molgraph_stats: boolean;
  calculate_molcomplexity: boolean;
  calculate_NIBR: boolean;
  nibr_scheduler: string;
  calculate_bredt: boolean;
  calculate_lilly: boolean;
  lilly_scheduler: string;
}

export interface SynthesisConfig {
  run: boolean;
  filter_solved_only: boolean;
  sa_score_min: number;
  sa_score_max: number | string;
  syba_score_min: number;
  syba_score_max: number | string;
  ra_score_min: number;
  ra_score_max: number;
}

export interface DockingToolConfig {
  bin: string;
  autobox_ligand?: string;
  autobox_add?: number;
  center_x?: number;
  center_y?: number;
  center_z?: number;
  size_x?: number;
  size_y?: number;
  size_z?: number;
  cpu?: number;
  seed?: number;
  exhaustiveness?: number;
  num_modes?: number;
  env_path?: string;
}

export interface DockingConfig {
  run: boolean;
  tools: string;
  receptor_pdb: string;
  auto_run: boolean;
  run_in_background: boolean;
  smina_config: DockingToolConfig;
  gnina_config: DockingToolConfig;
}

// Pipeline types
export type StageStatus = 'pending' | 'running' | 'completed' | 'error' | 'skipped';

export interface StageInfo {
  name: string;
  displayName: string;
  status: StageStatus;
  progress: number;
  message?: string;
  startTime?: Date;
  endTime?: Date;
}

export interface PipelineState {
  jobId: string | null;
  isRunning: boolean;
  stages: Record<string, StageInfo>;
  logs: LogEntry[];
  error: string | null;
}

export interface LogEntry {
  timestamp: Date;
  level: 'info' | 'warn' | 'error' | 'debug';
  message: string;
}

// RPC types
export interface RpcRequest {
  jsonrpc: '2.0';
  id: number;
  method: string;
  params?: Record<string, unknown>;
}

export interface RpcResponse {
  jsonrpc: '2.0';
  id: number;
  result?: unknown;
  error?: {
    code: number;
    message: string;
    data?: unknown;
  };
}

export interface RpcNotification {
  jsonrpc: '2.0';
  method: string;
  params: Record<string, unknown>;
}

// Validation
export interface ValidationResult {
  valid: boolean;
  errors: string[];
}

// File browser
export interface FileInfo {
  name: string;
  path: string;
  isDirectory: boolean;
  size?: number;
  modified?: Date;
}
