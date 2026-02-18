// Screen types
export type Screen =
  | 'welcome'
  | 'configMain'
  | 'configDescriptors'
  | 'configFilters'
  | 'configSynthesis'
  | 'configRetrosynthesis'
  | 'configDocking'
  | 'pipelineRunner'
  | 'history'
  | 'results'
  | 'wizardInputSelection'
  | 'wizardStageSelection'
  | 'wizardStageOrder'
  | 'wizardConfigMolPrep'
  | 'wizardConfigDescriptors'
  | 'wizardConfigFilters'
  | 'wizardConfigSynthesis'
  | 'wizardConfigDocking'
  | 'wizardConfigDockingFilters'
  | 'wizardReview';

// Screen shortcut definition
export interface ScreenShortcut {
  key: string;
  label: string;
  disabled?: boolean;
}

// Job status type
export type JobStatus = 'running' | 'completed' | 'cancelled' | 'error';

// Job history record
export interface JobHistoryRecord {
  id: string;
  name: string;
  startTime: string;
  endTime?: string;
  status: JobStatus;
  config: {
    inputPath: string;
    outputPath: string;
    stages: string[];
  };
  results?: {
    moleculesProcessed: number;
    moleculesFiltered?: number;
    dockingHits?: number;
  };
  error?: string;
}

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
  config_mol_prep: string;
  config_descriptors: string;
  config_structFilters: string;
  config_synthesis: string;
  config_docking: string;
  config_docking_filters: string;
}

export interface MolPrepConfig {
  run: boolean;
  n_jobs: number;
  steps: {
    keep_largest_fragment?: boolean;
    remove_stereochemistry?: boolean;
    fix_mol?: {
      enabled?: boolean;
    };
    remove_salts_solvents?: {
      enabled?: boolean;
    };
    standardize_mol?: {
      enabled?: boolean;
      uncharge?: boolean;
      stereo?: boolean;
    };
  };
  filters: {
    allowed_atoms?: string[];
    reject_radicals?: boolean;
    require_neutral?: boolean;
    reject_isotopes?: boolean;
    require_single_fragment?: boolean;
  };
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
  run_retrosynthesis: boolean;
  filter_solved_only: boolean;
  sa_score_min: number;
  sa_score_max: number | string;
  syba_score_min: number;
  syba_score_max: number | string;
  ra_score_min: number;
  ra_score_max: number;
}

export interface RetrosynthesisConfig {
  // AiZynthFinder expansion models
  expansion_uspto_model: string;
  expansion_uspto_templates: string;
  expansion_ringbreaker_model: string;
  expansion_ringbreaker_templates: string;
  // Filter model
  filter_uspto: string;
  // Stock database
  stock_zinc: string;
}

export interface DockingToolConfig {
  bin: string;
  env_path?: string;
  // Search space - autobox
  autobox_ligand?: string;
  autobox_add?: number;
  // Search space - manual
  center_x?: number;
  center_y?: number;
  center_z?: number;
  size_x?: number;
  size_y?: number;
  size_z?: number;
  // Flexible docking
  flex?: string;
  flexres?: string;
  flexdist_ligand?: string;
  flexdist?: number;
  // Scoring
  scoring?: string;
  custom_scoring?: string;
  custom_atoms?: string;
  // Modes
  score_only?: boolean;
  local_only?: boolean;
  minimize?: boolean;
  minimize_iters?: number;
  randomize_only?: boolean;
  // Output options
  energy_range?: number;
  min_rmsd_filter?: number;
  out_flex?: string;
  log?: string;
  quiet?: boolean;
  addH?: boolean;
  // Processing
  cpu?: number;
  seed?: number;
  exhaustiveness?: number;
  num_modes?: number;
  output_dir?: string;
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

export interface DockingFiltersConfig {
  run: boolean;
  run_after_docking?: boolean;
  input_sdf?: string | null;
  receptor_pdb?: string | null;
  pose_quality?: {
    enabled?: boolean;
    max_clashes?: number;
    max_strain_energy?: number;
    strain_forcefield?: string;
    clash_tolerance?: number;
  };
  interactions?: {
    enabled?: boolean;
    reference_ligand?: string | null;
    min_hbonds?: number;
    required_residues?: string[];
    forbidden_residues?: string[];
    interaction_types?: string[];
    similarity_threshold?: number;
  };
  shepherd_score?: {
    enabled?: boolean;
    reference_ligand?: string | null;
    min_shape_score?: number;
    alpha?: number;
    align_before_scoring?: boolean;
  };
  conformer_deviation?: {
    enabled?: boolean;
    num_conformers?: number;
    conformer_method?: string;
    max_rmsd_to_conformer?: number;
    random_seed?: number;
    optimize_conformers?: boolean;
  };
  aggregation?: {
    mode?: 'all' | 'any';
    save_metrics?: boolean;
    save_failed?: boolean;
  };
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

// Pipeline preflight
export type PreflightLevel = 'error' | 'warning' | 'info';

export interface PreflightCheck {
  code: string;
  level: PreflightLevel;
  stage?: string;
  field?: string;
  message: string;
}

export interface StagePreflightReport {
  stage: string;
  status: 'ok' | 'warning' | 'error';
  checks: PreflightCheck[];
}

export interface PipelinePreflightResult {
  valid: boolean;
  molecule_count: number | null;
  estimated_runtime: 'short' | 'medium' | 'long' | 'unknown';
  checks: PreflightCheck[];
  stage_reports: StagePreflightReport[];
}

// File browser
export interface FileInfo {
  name: string;
  path: string;
  isDirectory: boolean;
  size?: number;
  modified?: Date;
}

// Toast notification types
export type ToastType = 'success' | 'error' | 'info' | 'warning';

export interface Toast {
  id: string;
  type: ToastType;
  message: string;
  duration: number;
}

// Confirmation dialog types
export interface ConfirmDialogConfig {
  title: string;
  message: string;
  confirmLabel?: string;
  cancelLabel?: string;
  onConfirm: () => void;
  onCancel?: () => void;
}

// Help overlay types
export interface HelpShortcut {
  key: string;
  description: string;
}

export interface ScreenHelp {
  title: string;
  shortcuts: HelpShortcut[];
}

// Pipeline Wizard types
export type WizardStep =
  | 'stage-selection'
  | 'stage-order'
  | 'config-mol-prep'
  | 'config-descriptors'
  | 'config-filters'
  | 'config-synthesis'
  | 'config-docking'
  | 'review';

export interface WizardStageConfig {
  enabled: boolean;
  order: number;
  quickParams: Record<string, unknown>;
  preset?: string;
}

export interface WizardDependencies {
}

export interface WizardState {
  currentStep: WizardStep;
  selectedStages: string[];
  stageOrder: string[];
  stageConfigs: Record<string, WizardStageConfig>;
  preflight: PipelinePreflightResult | null;
  dependencies: WizardDependencies;
}

// Quick config parameter definitions
export interface QuickParamDef {
  key: string;
  label: string;
  type: 'boolean' | 'number' | 'range' | 'select';
  min?: number;
  max?: number;
  options?: string[];
  description?: string;
}

export interface StageQuickParams {
  stageName: string;
  displayName: string;
  params: QuickParamDef[];
  presets?: Record<string, Record<string, unknown>>;
}
