import React from 'react';
import Plot from 'react-plotly.js';
import { Data, Layout } from 'plotly.js';
import { DataPoint } from '../types';
import { PlotConfig } from './PlotControls';

interface PlotsProps {
  data: DataPoint[];
  plotConfig: PlotConfig;
}

const Plots: React.FC<PlotsProps> = ({ data, plotConfig }) => {
  const { scenarios, yAxis, xAxis } = plotConfig;

  // Filter data for selected scenarios
  const filteredData = data.filter(d => scenarios.includes(d.scenario));

  // Group data by scenario
  const dataByScenario = filteredData.reduce((acc, d) => {
    if (!acc[d.scenario]) {
      acc[d.scenario] = [];
    }
    acc[d.scenario].push(d);
    return acc;
  }, {} as Record<string, DataPoint[]>);

  // Get y-axis value based on selection
  const getYValue = (point: DataPoint) => {
    switch (yAxis) {
      case 'cost':
        return point.expectedCost;
      case 'expectedEffectiveness':
        return point.expectedEffectiveness;
      case 'effectiveCost':
        return point.effectiveCost;
      default:
        return 0;
    }
  };

  // Get x-axis value based on selection
  const getXValue = (point: DataPoint) => {
    return point[xAxis as keyof DataPoint] as number;
  };

  const plotData: { data: Data[]; layout: Partial<Layout> } = {
    data: Object.entries(dataByScenario).map(([scenario, points]) => ({
      x: points.map(getXValue),
      y: points.map(getYValue),
      type: 'scatter' as const,
      mode: 'lines' as const,
      name: scenario,
      line: { shape: 'spline' }
    })),
    layout: {
      title: `${yAxis} vs ${xAxis}`,
      xaxis: { title: xAxis },
      yaxis: { title: yAxis },
      font: {size: 22},
      showlegend: true
    }
  };

  return (
    <div>
      <Plot
        data={plotData.data}
        layout={plotData.layout}
        style={{ width: '100%', height: '400px' }}
      />
    </div>
  );
};

export default Plots; 
