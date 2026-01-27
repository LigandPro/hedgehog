import { useState, useCallback } from 'react';

export type FieldType = 'text' | 'number' | 'boolean' | 'select' | 'path';

export interface UseEditableFieldOptions {
  type: FieldType;
  initialValue: unknown;
  options?: string[]; // For select type
  onSave: (value: unknown) => void;
  onCancel?: () => void;
}

export interface UseEditableFieldReturn {
  isEditing: boolean;
  editValue: string;
  setEditValue: (value: string) => void;
  startEdit: () => void;
  cancelEdit: () => void;
  saveEdit: () => void;
  toggleBoolean: () => void;
  cycleSelect: () => void;
  handleInput: (input: string, key: {
    escape?: boolean;
    return?: boolean;
    backspace?: boolean;
    delete?: boolean;
  }) => boolean;
}

/**
 * Unified hook for editing field values
 *
 * Consistent keyboard handling:
 * - 'e' or Enter: Start edit mode
 * - Escape: Cancel edit
 * - Enter (in edit mode): Save value
 * - Space: Toggle boolean
 */
export function useEditableField({
  type,
  initialValue,
  options,
  onSave,
  onCancel,
}: UseEditableFieldOptions): UseEditableFieldReturn {
  const [isEditing, setIsEditing] = useState(false);
  const [editValue, setEditValue] = useState('');

  const startEdit = useCallback(() => {
    if (type === 'boolean') {
      // Toggle boolean directly
      onSave(!initialValue);
      return;
    }

    if (type === 'select' && options) {
      // Cycle select options
      const currentIdx = options.indexOf(String(initialValue));
      const nextIdx = (currentIdx + 1) % options.length;
      onSave(options[nextIdx]);
      return;
    }

    // Text/number: enter edit mode
    setEditValue(initialValue !== undefined && initialValue !== null ? String(initialValue) : '');
    setIsEditing(true);
  }, [type, initialValue, options, onSave]);

  const cancelEdit = useCallback(() => {
    setIsEditing(false);
    setEditValue('');
    onCancel?.();
  }, [onCancel]);

  const saveEdit = useCallback(() => {
    let finalValue: unknown = editValue;

    if (type === 'number') {
      const parsed = parseFloat(editValue);
      finalValue = isNaN(parsed) ? undefined : parsed;
    }

    onSave(finalValue);
    setIsEditing(false);
    setEditValue('');
  }, [type, editValue, onSave]);

  const toggleBoolean = useCallback(() => {
    if (type === 'boolean') {
      onSave(!initialValue);
    }
  }, [type, initialValue, onSave]);

  const cycleSelect = useCallback(() => {
    if (type === 'select' && options) {
      const currentIdx = options.indexOf(String(initialValue));
      const nextIdx = (currentIdx + 1) % options.length;
      onSave(options[nextIdx]);
    }
  }, [type, initialValue, options, onSave]);

  /**
   * Handle keyboard input for edit mode
   * Returns true if the input was handled (consumed)
   */
  const handleInput = useCallback((
    input: string,
    key: {
      escape?: boolean;
      return?: boolean;
      backspace?: boolean;
      delete?: boolean;
    }
  ): boolean => {
    if (!isEditing) {
      // Start edit on 'e' or Enter
      if (input === 'e' || key.return) {
        startEdit();
        return true;
      }
      // Toggle boolean on Space
      if (input === ' ' && type === 'boolean') {
        toggleBoolean();
        return true;
      }
      return false;
    }

    // In edit mode
    if (key.escape) {
      cancelEdit();
      return true;
    }

    if (key.return) {
      saveEdit();
      return true;
    }

    if (key.backspace || key.delete) {
      setEditValue((prev) => prev.slice(0, -1));
      return true;
    }

    // Add character
    if (input && input.length === 1) {
      // Number validation
      if (type === 'number') {
        if (!/[\d.\-]/.test(input)) {
          return true; // Consume but don't add invalid char
        }
      }
      setEditValue((prev) => prev + input);
      return true;
    }

    return true; // Consume all input in edit mode
  }, [isEditing, type, startEdit, cancelEdit, saveEdit, toggleBoolean]);

  return {
    isEditing,
    editValue,
    setEditValue,
    startEdit,
    cancelEdit,
    saveEdit,
    toggleBoolean,
    cycleSelect,
    handleInput,
  };
}

export default useEditableField;
