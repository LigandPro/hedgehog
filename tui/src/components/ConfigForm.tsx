import React, { useState } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { useTerminalSize } from '../hooks/useTerminalSize.js';

export interface FormField {
  key: string;
  label: string;
  type: 'text' | 'number' | 'boolean' | 'path';
  value: string | number | boolean;
  description?: string;
}

interface ConfigFormProps {
  fields: FormField[];
  onSave: (values: Record<string, string | number | boolean>) => void;
  onCancel: () => void;
}

export function ConfigForm({ fields, onSave, onCancel }: ConfigFormProps): React.ReactElement {
  const { width: terminalWidth } = useTerminalSize();

  const [selectedIndex, setSelectedIndex] = useState(0);
  const [editMode, setEditMode] = useState(false);
  const [values, setValues] = useState<Record<string, string | number | boolean>>(() => {
    const initial: Record<string, string | number | boolean> = {};
    for (const field of fields) {
      initial[field.key] = field.value;
    }
    return initial;
  });
  const [editValue, setEditValue] = useState('');

  useInput((input, key) => {
    if (editMode) {
      if (key.escape) {
        setEditMode(false);
      } else if (key.return) {
        const field = fields[selectedIndex];
        let newValue: string | number | boolean = editValue;
        
        if (field.type === 'number') {
          newValue = parseFloat(editValue) || 0;
        } else if (field.type === 'boolean') {
          newValue = editValue.toLowerCase() === 'true' || editValue === '1';
        }
        
        setValues({ ...values, [field.key]: newValue });
        setEditMode(false);
      }
      return;
    }

    if (key.upArrow) {
      setSelectedIndex(Math.max(0, selectedIndex - 1));
    } else if (key.downArrow) {
      setSelectedIndex(Math.min(fields.length - 1, selectedIndex + 1));
    } else if (key.return || input === 'e') {
      const field = fields[selectedIndex];
      if (field.type === 'boolean') {
        setValues({ ...values, [field.key]: !values[field.key] });
      } else {
        setEditValue(String(values[field.key]));
        setEditMode(true);
      }
    } else if (input === 's') {
      onSave(values);
    } else if (key.escape || input === 'q') {
      onCancel();
    }
  });

  return (
    <Box flexDirection="column">
      {fields.map((field, index) => {
        const isSelected = index === selectedIndex;
        const isEditing = isSelected && editMode;
        const value = values[field.key];

        return (
          <Box key={field.key} marginY={0}>
            <Text color={isSelected ? 'cyan' : 'white'}>
              {isSelected ? '▶ ' : '  '}
            </Text>
            <Text dimColor>{field.label.padEnd(25)}</Text>
            {isEditing ? (
              <Box>
                <TextInput
                  value={editValue}
                  onChange={setEditValue}
                  focus={true}
                />
              </Box>
            ) : (
              <Text color={field.type === 'boolean' ? (value ? 'green' : 'red') : 'yellow'}>
                {field.type === 'boolean' ? (value ? 'Yes' : 'No') : String(value)}
              </Text>
            )}
          </Box>
        );
      })}
      <Box marginTop={1}>
        <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
      </Box>
      <Box gap={2}>
        <Text><Text color="cyan">[s]</Text> Save</Text>
        <Text><Text color="cyan">[e/Enter]</Text> Edit</Text>
        <Text><Text color="cyan">[Esc]</Text> Cancel</Text>
      </Box>
    </Box>
  );
}

export default ConfigForm;
