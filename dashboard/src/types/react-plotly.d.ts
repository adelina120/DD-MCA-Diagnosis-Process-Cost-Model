declare module 'react-plotly.js' {
  import { Component } from 'react';
  import { Data, Layout, Config } from 'plotly.js';

  interface PlotParams {
    data: Data[];
    layout?: Partial<Layout>;
    config?: Partial<Config>;
    style?: React.CSSProperties;
    className?: string;
    onInitialized?: (figure: any, graphDiv: HTMLElement) => void;
    onUpdate?: (figure: any, graphDiv: HTMLElement) => void;
    onPurge?: (figure: any, graphDiv: HTMLElement) => void;
    onError?: (err: Error) => void;
    onAfterPlot?: () => void;
    onSelected?: (event: any) => void;
    onSelecting?: (event: any) => void;
    onUnselect?: () => void;
    onRelayout?: (event: any) => void;
    onRedraw?: () => void;
    onAnimated?: () => void;
    onAnimatingFrame?: (event: any) => void;
    onDoubleClick?: () => void;
    onHover?: (event: any) => void;
    onUnhover?: (event: any) => void;
    onClick?: (event: any) => void;
    onDeselect?: () => void;
    onRelayouting?: (event: any) => void;
    onRestyle?: (event: any) => void;
    onRestyling?: (event: any) => void;
    onSliderChange?: (event: any) => void;
    onSliderEnd?: (event: any) => void;
    onSliderStart?: (event: any) => void;
    onTransitioning?: (event: any) => void;
    onTransitionInterrupted?: (event: any) => void;
    onSunburstClick?: (event: any) => void;
    onSunburstHover?: (event: any) => void;
    onSunburstUnhover?: (event: any) => void;
    onSunburstSelect?: (event: any) => void;
    onSunburstDeselect?: (event: any) => void;
    onSunburstRestyle?: (event: any) => void;
    onSunburstRelayout?: (event: any) => void;
    onSunburstDoubleClick?: (event: any) => void;
    onSunburstTransitioning?: (event: any) => void;
    onSunburstTransitionInterrupted?: (event: any) => void;
    onSunburstSliderChange?: (event: any) => void;
    onSunburstSliderEnd?: (event: any) => void;
    onSunburstSliderStart?: (event: any) => void;
    onSunburstAnimatingFrame?: (event: any) => void;
    onSunburstAnimated?: () => void;
    onSunburstRedraw?: () => void;
    onSunburstRestyling?: (event: any) => void;
    onSunburstRelayouting?: (event: any) => void;
    onSunburstSelecting?: (event: any) => void;
    onSunburstUnselect?: () => void;
    onSunburstHovering?: (event: any) => void;
    onSunburstUnhovering?: (event: any) => void;
    onSunburstClicking?: (event: any) => void;
    onSunburstDoubleClicking?: (event: any) => void;
  }

  export default class Plot extends Component<PlotParams> {}
} 