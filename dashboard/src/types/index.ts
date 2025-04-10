export interface Parameters {
  n: number;
  cmaYield: number;
  gpYield: number;
  wesYield1Tier: number;
  wesYield2Tier: number;
  wesYield3Tier: number;
  expertFee: number;
  cmaCost: number;
  gpCost: number;
  wesCost: number;
  cmaPPV: number;
  cmaNPV: number;
  gpPPV: number;
  gpNPV: number;
  wesPPV: number;
  wesNPV: number;
  aiPrecision: number;
  aiFDR: number;
  aiFOR: number;
  aiNPV: number;
  alpha: number;
  lambda: number;
}

export interface DataPoint {
  scenario: string;
  expectedCost: number;
  expectedUtility: number;
  effectiveCost: number;
  cmaCost: number;
  gpCost: number;
  aiPerformance: number;
  alphaValues: number;
}

export interface PlotData {
  data: DataPoint[];
  layout: any;
  config: any;
} 