import React from 'react';
import { Box, FormControl, InputLabel, Select, MenuItem, Typography, SelectChangeEvent } from '@mui/material';
import { Parameters } from '../types';

interface PlotControlsProps {
  parameters: Parameters;
  onPlotConfigChange: (config: PlotConfig) => void;
}

export interface PlotConfig {
  scenarios: string[];
  yAxis: 'cost' | 'expectedQALY' | 'costPerQALY';
  xAxis: keyof Parameters;
  showAIComparison: boolean;
}

const availableScenarios = [
  "Expert-alone: Scenario 1 (CMA + GP)",
  "Expert-alone: Scenario 2 (CMA + GP + WES)",
  "Expert-alone: Scenario 3 (CMA + WES)",
  "Expert-alone: Scenario 4 (WES alone)",
  "AI-delegation: r>r* (CMA + GP + WES)"
];

const yAxisOptions = [
  { value: 'cost', label: 'Cost' },
  { value: 'expectedQALY', label: 'Expected QALY' },
  { value: 'costPerQALY', label: 'Cost Per QALY' }
];

const PlotControls: React.FC<PlotControlsProps> = ({ parameters, onPlotConfigChange }) => {
  const [config, setConfig] = React.useState<PlotConfig>({
    scenarios: availableScenarios,
    yAxis: 'costPerQALY',
    xAxis: 'numberOfYears',
    showAIComparison: false
  });

  const handleScenarioChange = (event: SelectChangeEvent<string[]>) => {
    const newConfig = {
      ...config,
      scenarios: event.target.value as string[]
    };
    setConfig(newConfig);
    onPlotConfigChange(newConfig);
  };

  const handleYAxisChange = (event: SelectChangeEvent<'cost' | 'expectedQALY' | 'costPerQALY'>) => {
    const newConfig = {
      ...config,
      yAxis: event.target.value as 'cost' | 'expectedQALY' | 'costPerQALY'
    };
    setConfig(newConfig);
    onPlotConfigChange(newConfig);
  };

  const handleXAxisChange = (event: SelectChangeEvent<keyof Parameters>) => {
    const newConfig = {
      ...config,
      xAxis: event.target.value as keyof Parameters
    };
    setConfig(newConfig);
    onPlotConfigChange(newConfig);
  };

  // Get available x-axis parameters based on selected scenarios and y-axis
  const getAvailableXAxisParams = () => {
    const baseParams = [
      'cmaYield',
      'gpYield',
      'wesYield1Tier',
      'wesYield2Tier',
      'wesYield3Tier',
      'cmaPPV',
      'cmaNPV',
      'gpPPV',
      'gpNPV',
      'wesPPV',
      'wesNPV',
      'cmaTAT',
      'gpTAT',
      'wesTAT',
      'expertFee',
      'expertTAT',
      'cmaCost',
      'gpCost',
      'wesCost',
    ];

    // Add AI-related parameters if AI scenario is selected
    if (config.scenarios.includes("AI-delegation: r>r* (CMA + GP + WES)")) {
      baseParams.push('aiPrecision', 'aiFDR', 'aiFOR', 'aiNPV');
    }

    // Add utility measurements if y-axis is not cost
    if (config.yAxis !== 'cost') {
      baseParams.push('uTP','uFP','uTN','uFN','uInitial','numberOfYears');
    }

    return baseParams;
  };

  return (
    <Box sx={{ mb: 2 }}>
      <Typography gutterBottom>Plot Configuration</Typography>
      
      <FormControl fullWidth margin="normal">
        <InputLabel>Scenarios</InputLabel>
        <Select
          multiple
          value={config.scenarios}
          onChange={handleScenarioChange}
          label="Scenarios"
        >
          {availableScenarios.map(scenario => (
            <MenuItem key={scenario} value={scenario}>
              {scenario}
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <FormControl fullWidth margin="normal">
        <InputLabel>Y-Axis</InputLabel>
        <Select
          value={config.yAxis}
          onChange={handleYAxisChange}
          label="Y-Axis"
        >
          {yAxisOptions.map(option => (
            <MenuItem key={option.value} value={option.value}>
              {option.label}
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <FormControl fullWidth margin="normal">
        <InputLabel>X-Axis</InputLabel>
        <Select
          value={config.xAxis}
          onChange={handleXAxisChange}
          label="X-Axis"
        >
          {getAvailableXAxisParams().map(param => (
            <MenuItem key={param} value={param}>
              {param}
            </MenuItem>
          ))}
        </Select>
      </FormControl>
    </Box>
  );
};

export default PlotControls; 
