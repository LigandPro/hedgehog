import { useState, useCallback, useMemo } from 'react';

export interface UseSearchOptions<T> {
  items: T[];
  searchFields: (item: T) => string[];
  enabled?: boolean;
}

export interface UseSearchReturn<T> {
  searchQuery: string;
  setSearchQuery: (query: string) => void;
  searchActive: boolean;
  setSearchActive: (active: boolean) => void;
  filteredItems: T[];
  filteredIndices: number[];
  highlightMatch: (text: string) => { text: string; highlighted: boolean }[];
  clearSearch: () => void;
  handleSearchInput: (input: string, key: { escape?: boolean; return?: boolean }) => boolean;
}

/**
 * Hook for search/filter functionality in list screens
 *
 * Usage:
 * 1. Press '/' to activate search mode
 * 2. Type to filter items
 * 3. Press Escape to clear and exit search mode
 * 4. Press Enter to confirm and exit search mode (keeping filter)
 */
export function useSearch<T>({
  items,
  searchFields,
  enabled = true,
}: UseSearchOptions<T>): UseSearchReturn<T> {
  const [searchQuery, setSearchQuery] = useState('');
  const [searchActive, setSearchActive] = useState(false);

  const clearSearch = useCallback(() => {
    setSearchQuery('');
    setSearchActive(false);
  }, []);

  // Fuzzy match function
  const fuzzyMatch = useCallback((text: string, query: string): boolean => {
    if (!query) return true;
    const lowerText = text.toLowerCase();
    const lowerQuery = query.toLowerCase();

    // Simple substring match
    if (lowerText.includes(lowerQuery)) return true;

    // Fuzzy match: all query chars must appear in order
    let queryIdx = 0;
    for (let i = 0; i < lowerText.length && queryIdx < lowerQuery.length; i++) {
      if (lowerText[i] === lowerQuery[queryIdx]) {
        queryIdx++;
      }
    }
    return queryIdx === lowerQuery.length;
  }, []);

  // Filter items based on search query
  const { filteredItems, filteredIndices } = useMemo(() => {
    if (!searchQuery || !enabled) {
      return {
        filteredItems: items,
        filteredIndices: items.map((_, i) => i),
      };
    }

    const filtered: T[] = [];
    const indices: number[] = [];

    items.forEach((item, index) => {
      const fields = searchFields(item);
      const matches = fields.some((field) => fuzzyMatch(field, searchQuery));
      if (matches) {
        filtered.push(item);
        indices.push(index);
      }
    });

    return { filteredItems: filtered, filteredIndices: indices };
  }, [items, searchQuery, searchFields, fuzzyMatch, enabled]);

  // Highlight matching parts of text
  const highlightMatch = useCallback((text: string): { text: string; highlighted: boolean }[] => {
    if (!searchQuery) {
      return [{ text, highlighted: false }];
    }

    const lowerText = text.toLowerCase();
    const lowerQuery = searchQuery.toLowerCase();
    const matchIndex = lowerText.indexOf(lowerQuery);

    if (matchIndex === -1) {
      return [{ text, highlighted: false }];
    }

    const parts: { text: string; highlighted: boolean }[] = [];

    if (matchIndex > 0) {
      parts.push({ text: text.slice(0, matchIndex), highlighted: false });
    }

    parts.push({
      text: text.slice(matchIndex, matchIndex + searchQuery.length),
      highlighted: true,
    });

    if (matchIndex + searchQuery.length < text.length) {
      parts.push({
        text: text.slice(matchIndex + searchQuery.length),
        highlighted: false,
      });
    }

    return parts;
  }, [searchQuery]);

  /**
   * Handle keyboard input for search mode
   * Returns true if the input was handled (consumed)
   */
  const handleSearchInput = useCallback((
    input: string,
    key: { escape?: boolean; return?: boolean; backspace?: boolean; delete?: boolean }
  ): boolean => {
    if (!enabled) return false;

    // Activate search mode with '/'
    if (!searchActive && input === '/') {
      setSearchActive(true);
      return true;
    }

    if (!searchActive) return false;

    // Exit search mode
    if (key.escape) {
      clearSearch();
      return true;
    }

    // Confirm search (exit search mode but keep filter)
    if (key.return) {
      setSearchActive(false);
      return true;
    }

    // Backspace
    if (key.backspace || key.delete) {
      setSearchQuery((prev) => prev.slice(0, -1));
      return true;
    }

    // Add character to search
    if (input && input.length === 1 && input !== '/') {
      setSearchQuery((prev) => prev + input);
      return true;
    }

    return true; // Consume all input in search mode
  }, [enabled, searchActive, clearSearch]);

  return {
    searchQuery,
    setSearchQuery,
    searchActive,
    setSearchActive,
    filteredItems,
    filteredIndices,
    highlightMatch,
    clearSearch,
    handleSearchInput,
  };
}

export default useSearch;
