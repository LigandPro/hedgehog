export function getByPath<T = unknown>(obj: Record<string, unknown>, path?: string[]): T | undefined {
  if (!path || path.length === 0) return undefined;
  let val: unknown = obj;
  for (const key of path) {
    if (val === null || val === undefined) return undefined;
    val = (val as Record<string, unknown>)[key];
  }
  return val as T;
}

export function setByPath(
  obj: Record<string, unknown>,
  path: string[],
  value: unknown
): Record<string, unknown> {
  if (path.length === 0) return obj;

  const result = { ...obj };
  let current: Record<string, unknown> = result;

  for (let i = 0; i < path.length - 1; i++) {
    const key = path[i];
    current[key] = { ...(current[key] as Record<string, unknown> || {}) };
    current = current[key] as Record<string, unknown>;
  }

  current[path[path.length - 1]] = value;
  return result;
}
