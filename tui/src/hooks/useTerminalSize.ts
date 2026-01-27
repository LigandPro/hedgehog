import { useState, useEffect } from 'react';
import { useStdout } from 'ink';

interface TerminalSize {
  width: number;
  height: number;
}

/**
 * Hook that returns terminal dimensions and updates on resize.
 * Unlike useStdout alone, this triggers re-renders when terminal size changes.
 */
export function useTerminalSize(): TerminalSize {
  const { stdout } = useStdout();

  const [size, setSize] = useState<TerminalSize>({
    width: stdout?.columns || 80,
    height: stdout?.rows || 24,
  });

  useEffect(() => {
    if (!stdout) return;

    const handleResize = () => {
      setSize({
        width: stdout.columns || 80,
        height: stdout.rows || 24,
      });
    };

    // Listen for resize events
    stdout.on('resize', handleResize);

    // Initial size
    handleResize();

    return () => {
      stdout.off('resize', handleResize);
    };
  }, [stdout]);

  return size;
}

export default useTerminalSize;
