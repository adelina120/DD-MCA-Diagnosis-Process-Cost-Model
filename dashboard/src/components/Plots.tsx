import React from 'react';
import Plot from 'react-plotly.js';
import { Box, Typography } from '@mui/material';
import { DataPoint } from '../types';
import { Data, Layout } from 'plotly.js';

interface PlotsProps {
  data: DataPoint[];
}

const Plots: React.FC<PlotsProps> = ({ data }) => {
  const expertAloneData = data.filter(d => d.scenario !== "AI-delegation (r>r*)");
  const aiDelegationData = data.filter(d => 
    d.scenario === "Scenario 2 (CMA + GP + WES)" || 
    d.scenario === "Scenario 3 (CMA + WES)" || 
    d.scenario === "AI-delegation (r>r*)"
  );

  // Group data by scenario
  const expertAloneByScenario = expertAloneData.reduce((acc, d) => {
    if (!acc[d.scenario]) {
      acc[d.scenario] = [];
    }
    acc[d.scenario].push(d);
    return acc;
  }, {} as Record<string, DataPoint[]>);

  const aiDelegationByScenario = aiDelegationData.reduce((acc, d) => {
    if (!acc[d.scenario]) {
      acc[d.scenario] = [];
    }
    acc[d.scenario].push(d);
    return acc;
  }, {} as Record<string, DataPoint[]>);

  const effectiveCostPlot: { data: Data[]; layout: Partial<Layout> } = {
    data: Object.entries(expertAloneByScenario).map(([scenario, points]) => ({
      x: points.map(p => p.alphaValues),
      y: points.map(p => p.effectiveCost),
      type: 'scatter' as const,
      mode: 'lines' as const,
      name: scenario,
      line: { shape: 'spline' }
    })),
    layout: {
      title: 'Effective Cost vs Alpha',
      xaxis: { title: 'Alpha' },
      yaxis: { title: 'Effective Cost' },
      showlegend: true
    }
  };

  const expectedUtilityPlot: { data: Data[]; layout: Partial<Layout> } = {
    data: Object.entries(expertAloneByScenario).map(([scenario, points]) => ({
      x: points.map(p => p.alphaValues),
      y: points.map(p => p.expectedUtility),
      type: 'scatter' as const,
      mode: 'lines' as const,
      name: scenario,
      line: { shape: 'spline' }
    })),
    layout: {
      title: 'Expected Utility vs Alpha',
      xaxis: { title: 'Alpha' },
      yaxis: { title: 'Expected Utility' },
      showlegend: true
    }
  };

  const aiComparisonPlot: { data: Data[]; layout: Partial<Layout> } = {
    data: Object.entries(aiDelegationByScenario).map(([scenario, points]) => ({
      x: points.map(p => p.alphaValues),
      y: points.map(p => p.expectedUtility),
      type: 'scatter' as const,
      mode: 'lines' as const,
      name: scenario,
      line: { shape: 'spline' }
    })),
    layout: {
      title: 'AI vs Expert Comparison',
      xaxis: { title: 'Alpha' },
      yaxis: { title: 'Expected Utility' },
      showlegend: true
    }
  };

  return (
    <Box>
      <Typography variant="h6" gutterBottom>
        Model Visualizations
      </Typography>
      <Box sx={{ mb: 4 }}>
        <Plot
          data={effectiveCostPlot.data}
          layout={effectiveCostPlot.layout}
          style={{ width: '100%', height: '400px' }}
        />
      </Box>
      <Box sx={{ mb: 4 }}>
        <Plot
          data={expectedUtilityPlot.data}
          layout={expectedUtilityPlot.layout}
          style={{ width: '100%', height: '400px' }}
        />
      </Box>
      <Box>
        <Plot
          data={aiComparisonPlot.data}
          layout={aiComparisonPlot.layout}
          style={{ width: '100%', height: '400px' }}
        />
      </Box>
    </Box>
  );
};

export default Plots; 