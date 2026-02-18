import type { Screen } from '../../types/index.js';

export type WizardStageName =
  | 'mol_prep'
  | 'descriptors'
  | 'struct_filters'
  | 'synthesis'
  | 'docking'
  | 'docking_filters';

export interface StageMetadata {
  title: string;
  shortDescription: string;
  whatItDoes: string;
  reads: string[];
  writes: string[];
  heavyLevel: 'low' | 'medium' | 'high';
  keyParamsMap: Record<string, string>;
  configScreen: Screen;
}

export const WIZARD_STAGE_ORDER: WizardStageName[] = [
  'mol_prep',
  'descriptors',
  'struct_filters',
  'synthesis',
  'docking',
  'docking_filters',
];

export const STAGE_METADATA: Record<WizardStageName, StageMetadata> = {
  mol_prep: {
    title: 'Mol Prep',
    shortDescription: 'Standardize molecules (Datamol)',
    whatItDoes: 'Normalizes and filters molecules before downstream stages.',
    reads: ['Input molecules from main config'],
    writes: ['Standardized molecule set for next stages'],
    heavyLevel: 'medium',
    keyParamsMap: {
      n_jobs: 'Parallel jobs',
      require_neutral: 'Require neutral',
      require_single_fragment: 'Single fragment only',
    },
    configScreen: 'wizardConfigMolPrep',
  },
  descriptors: {
    title: 'Descriptors',
    shortDescription: 'Calculate molecular descriptors',
    whatItDoes: 'Calculates physicochemical descriptors and can filter by borders.',
    reads: ['Prepared molecules'],
    writes: ['Descriptor columns and filtered subset'],
    heavyLevel: 'medium',
    keyParamsMap: {
      batch_size: 'Batch size',
      filter_data: 'Filter output',
      molWt_min: 'MolWt min',
      molWt_max: 'MolWt max',
      logP_min: 'logP min',
      logP_max: 'logP max',
    },
    configScreen: 'wizardConfigDescriptors',
  },
  struct_filters: {
    title: 'Struct Filters',
    shortDescription: 'Apply structural filters',
    whatItDoes: 'Applies medicinal-chemistry alerts and structural filters.',
    reads: ['Prepared/descriptor data'],
    writes: ['Flags and optionally filtered molecules'],
    heavyLevel: 'medium',
    keyParamsMap: {
      calculate_common_alerts: 'Common alerts',
      calculate_NIBR: 'NIBR filters',
      calculate_lilly: 'Lilly filters',
      filter_data: 'Filter output',
    },
    configScreen: 'wizardConfigFilters',
  },
  synthesis: {
    title: 'Synthesis',
    shortDescription: 'Score synthesizability',
    whatItDoes: 'Computes SA/SYBA/RA and optionally runs retrosynthesis.',
    reads: ['Filtered molecules'],
    writes: ['Synthesis scores and optional route data'],
    heavyLevel: 'high',
    keyParamsMap: {
      run_retrosynthesis: 'Run retrosynthesis',
      sa_score_min: 'SA min',
      sa_score_max: 'SA max',
      ra_score_min: 'RA min',
      ra_score_max: 'RA max',
    },
    configScreen: 'wizardConfigSynthesis',
  },
  docking: {
    title: 'Docking',
    shortDescription: 'Run molecular docking',
    whatItDoes: 'Generates binding poses and docking scores against the receptor.',
    reads: ['Ligands and receptor PDB'],
    writes: ['Docking poses and score tables'],
    heavyLevel: 'high',
    keyParamsMap: {
      tools: 'Tool',
      exhaustiveness: 'Exhaustiveness',
      num_modes: 'Num modes',
    },
    configScreen: 'wizardConfigDocking',
  },
  docking_filters: {
    title: 'Docking Filters',
    shortDescription: 'Filter docking poses and interactions',
    whatItDoes: 'Filters docking output by geometry, interactions, and conformer checks.',
    reads: ['Docking poses and receptor context'],
    writes: ['Filtered docking candidates and metrics'],
    heavyLevel: 'high',
    keyParamsMap: {
      aggregation_mode: 'Aggregation mode',
      max_clashes: 'Max clashes',
      min_hbonds: 'Min H-bonds',
      max_rmsd_to_conformer: 'Max RMSD',
    },
    configScreen: 'wizardConfigDockingFilters',
  },
};

export function getStageMetadata(stage: string): StageMetadata | undefined {
  return STAGE_METADATA[stage as WizardStageName];
}

export function getStageScreen(stage: string): Screen {
  return getStageMetadata(stage)?.configScreen ?? 'wizardReview';
}

function formatValue(value: unknown): string {
  if (typeof value === 'boolean') return value ? 'Yes' : 'No';
  if (value === undefined || value === null || value === '') return '-';
  return String(value);
}

export function getStageSummary(
  stageName: string,
  quickParams: Record<string, unknown>,
  preset?: string
): string {
  switch (stageName) {
    case 'mol_prep':
      return 'Datamol standardization';
    case 'descriptors': {
      const batch = formatValue(quickParams.batch_size);
      return preset ? `Batch: ${batch} | Preset: ${preset}` : `Batch: ${batch}`;
    }
    case 'struct_filters':
      return `NIBR: ${formatValue(quickParams.calculate_NIBR)} | Lilly: ${formatValue(quickParams.calculate_lilly)}`;
    case 'synthesis':
      return `SA: ${formatValue(quickParams.sa_score_min)}-${formatValue(quickParams.sa_score_max)} | RA: ${formatValue(quickParams.ra_score_min)}-${formatValue(quickParams.ra_score_max)}`;
    case 'docking':
      return `Tool: ${formatValue(quickParams.tools)} | Exhaust: ${formatValue(quickParams.exhaustiveness)} | Modes: ${formatValue(quickParams.num_modes)}`;
    case 'docking_filters':
      return `Mode: ${formatValue(quickParams.aggregation_mode)} | Clashes<=${formatValue(quickParams.max_clashes)} | H-bonds>=${formatValue(quickParams.min_hbonds)} | RMSD<=${formatValue(quickParams.max_rmsd_to_conformer)}`;
    default:
      return 'Configured';
  }
}

export function getStageKeyParams(
  stageName: string,
  quickParams: Record<string, unknown>
): string[] {
  const metadata = getStageMetadata(stageName);
  if (!metadata) return [];

  return Object.entries(metadata.keyParamsMap)
    .filter(([param]) => quickParams[param] !== undefined)
    .map(([param, label]) => `${label}: ${formatValue(quickParams[param])}`);
}
