import { useEffect } from 'react';
import { useStore } from '../store/index.js';

export function useInputLock(active: boolean): void {
  const setInputLocked = useStore((state) => state.setInputLocked);

  useEffect(() => {
    if (!active) return;
    setInputLocked(true);
    return () => {
      setInputLocked(false);
    };
  }, [active, setInputLocked]);
}

export default useInputLock;
