import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { homedir } from 'node:os';
import path from 'node:path';
import React, { createContext, useCallback, useContext, useMemo, useState } from 'react';
import { DEFAULT_THEME_ID, THEMES } from './themes.js';
import type { ThemeDefinition, ThemeId } from './themes.js';

interface ThemeContextValue {
  theme: ThemeDefinition;
  themeId: ThemeId;
  themeIds: ThemeId[];
  setThemeId: (next: ThemeId) => void;
  setPreviewThemeId: (next: ThemeId | null) => void;
  clearPreviewTheme: () => void;
  cycleTheme: () => void;
}

const ThemeContext = createContext<ThemeContextValue | undefined>(undefined);
const THEME_CONFIG_PATH = path.join(homedir(), '.config', 'hedgehog', 'tui-theme.json');

interface ThemeConfig {
  themeId?: string;
}

const LEGACY_THEME_ALIASES = {
  opencode: 'ligadpro',
} as const;

function isThemeId(value: string): value is ThemeId {
  return Object.prototype.hasOwnProperty.call(THEMES, value);
}

function resolveThemeId(value: string): ThemeId | null {
  if (isThemeId(value)) return value;
  const aliased = LEGACY_THEME_ALIASES[value as keyof typeof LEGACY_THEME_ALIASES];
  if (aliased && isThemeId(aliased)) {
    return aliased;
  }
  return null;
}

function loadThemeId(): ThemeId {
  try {
    const content = readFileSync(THEME_CONFIG_PATH, 'utf8');
    const parsed = JSON.parse(content) as ThemeConfig;
    if (typeof parsed.themeId === 'string') {
      const resolved = resolveThemeId(parsed.themeId);
      if (resolved) {
        return resolved;
      }
    }
  } catch {
    // Fall back to the default theme when file is missing or invalid.
  }
  return DEFAULT_THEME_ID;
}

function saveThemeId(themeId: ThemeId): void {
  try {
    mkdirSync(path.dirname(THEME_CONFIG_PATH), { recursive: true });
    writeFileSync(THEME_CONFIG_PATH, JSON.stringify({ themeId }, null, 2), 'utf8');
  } catch {
    // Keep UI responsive even if persistence is unavailable.
  }
}

export function ThemeProvider({ children }: { children: React.ReactNode }): React.ReactElement {
  const [themeId, setThemeIdState] = useState<ThemeId>(() => loadThemeId());
  const [previewThemeId, setPreviewThemeIdState] = useState<ThemeId | null>(null);
  const themeIds = useMemo(() => Object.keys(THEMES) as ThemeId[], []);
  const persistedThemeId = THEMES[themeId] ? themeId : DEFAULT_THEME_ID;
  const activeThemeId = previewThemeId && THEMES[previewThemeId]
    ? previewThemeId
    : persistedThemeId;

  const setThemeId = useCallback((next: ThemeId) => {
    if (!THEMES[next]) return;
    setThemeIdState(next);
    setPreviewThemeIdState(null);
    saveThemeId(next);
  }, []);

  const setPreviewThemeId = useCallback((next: ThemeId | null) => {
    if (next === null) {
      setPreviewThemeIdState(null);
      return;
    }
    if (!THEMES[next]) return;
    setPreviewThemeIdState(next);
  }, []);

  const clearPreviewTheme = useCallback(() => {
    setPreviewThemeIdState(null);
  }, []);

  const cycleTheme = useCallback(() => {
    const currentIndex = themeIds.indexOf(persistedThemeId);
    const nextIndex = currentIndex === -1 ? 0 : (currentIndex + 1) % themeIds.length;
    const nextThemeId = themeIds[nextIndex] || DEFAULT_THEME_ID;
    setThemeIdState(nextThemeId);
    setPreviewThemeIdState(null);
    saveThemeId(nextThemeId);
  }, [persistedThemeId, themeIds]);

  const value = useMemo<ThemeContextValue>(() => ({
    theme: THEMES[activeThemeId],
    themeId: persistedThemeId,
    themeIds,
    setThemeId,
    setPreviewThemeId,
    clearPreviewTheme,
    cycleTheme,
  }), [
    activeThemeId,
    clearPreviewTheme,
    cycleTheme,
    persistedThemeId,
    setPreviewThemeId,
    setThemeId,
    themeIds,
  ]);

  return <ThemeContext.Provider value={value}>{children}</ThemeContext.Provider>;
}

export function useTheme(): ThemeContextValue {
  const value = useContext(ThemeContext);
  if (!value) {
    throw new Error('useTheme must be used inside ThemeProvider');
  }
  return value;
}
