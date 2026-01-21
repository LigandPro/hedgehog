import React, { useEffect } from 'react';
import { Box, useApp } from 'ink';
import { useStore } from './store/index.js';
import { getBridge, destroyBridge } from './services/python-bridge.js';
import { logger } from './utils/logger.js';

// Screens
import { Welcome } from './screens/Welcome.js';
import { ConfigMain } from './screens/ConfigMain.js';
import { ConfigDescriptors } from './screens/ConfigDescriptors.js';
import { ConfigFilters } from './screens/ConfigFilters.js';
import { ConfigSynthesis } from './screens/ConfigSynthesis.js';
import { ConfigDocking } from './screens/ConfigDocking.js';
import { PipelineRunner } from './screens/PipelineRunner.js';

export function App(): React.ReactElement {
  const { exit } = useApp();
  const screen = useStore((state) => state.screen);
  const setBackendReady = useStore((state) => state.setBackendReady);

  // Initialize Python backend
  useEffect(() => {
    const initBackend = async () => {
      try {
        const bridge = getBridge();
        await bridge.start();
        setBackendReady(true);
        logger.info('Backend initialized');
      } catch (error) {
        logger.error('Failed to initialize backend:', error);
        // Continue without backend - some features will be limited
      }
    };

    initBackend();

    return () => {
      destroyBridge();
    };
  }, []);

  // Render current screen
  const renderScreen = () => {
    switch (screen) {
      case 'welcome':
        return <Welcome />;
      case 'configMain':
        return <ConfigMain />;
      case 'configDescriptors':
        return <ConfigDescriptors />;
      case 'configFilters':
        return <ConfigFilters />;
      case 'configSynthesis':
        return <ConfigSynthesis />;
      case 'configDocking':
        return <ConfigDocking />;
      case 'pipelineRunner':
        return <PipelineRunner />;
      default:
        return <Welcome />;
    }
  };

  return (
    <Box flexDirection="column">
      {renderScreen()}
    </Box>
  );
}

export default App;
