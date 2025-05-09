# MCA Diagnosis Process Cost Model Dashboard

This dashboard provides an interactive interface for exploring and visualizing the MCA Diagnosis Process Cost Model. It allows users to adjust various parameters and see their impact on different scenarios in real-time.

## Features

- Interactive parameter controls for all model variables
- Real-time visualization of model outputs
- Comparison of different testing scenarios
- AI vs Expert performance comparison
- Responsive design for different screen sizes

## Prerequisites

- Node.js (version 14 or higher)
- npm (version 6 or higher)

## Installation

1. Navigate to the dashboard directory:
   ```bash
   cd dashboard
   ```

2. Install dependencies:
   ```bash
   npm install
   ```

## Running the Dashboard

1. Start the development server:
   ```bash
   npm start
   ```

2. Open your browser and navigate to:
   ```
   http://localhost:3000
   ```

## Building for Production

To create a production build:

```bash
npm run build
```

The build files will be created in the `build` directory.

## Project Structure

- `src/components/` - React components
  - `ParameterControls.tsx` - Controls for adjusting model parameters
  - `Plots.tsx` - Visualization components
- `src/types/` - TypeScript type definitions
- `src/utils/` - Utility functions
- `public/` - Static assets

## Available Parameters

The dashboard allows you to adjust the following parameters:

### Cost Parameters
- CMA Cost
- GP Cost
- WES Cost
- Expert Fee

### Test Performance Parameters
- CMA Yield
- CMA PPV
- CMA NPV
- GP Yield
- GP PPV
- GP NPV
- WES Yield (1st Tier)
- WES Yield (2nd Tier)
- WES Yield (3rd Tier)
- WES PPV
- WES NPV

### Utility Parameters
- Alpha
- Lambda

### AI Parameters
- AI Precision

## License

This project is licensed under the MIT License. 
