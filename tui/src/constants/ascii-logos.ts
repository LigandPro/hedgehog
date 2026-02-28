export interface AsciiLogo {
  logo: string;
}

// Progressive ASCII logos - from full HEDGEHOG down to compact fallback.
export const ASCII_LOGOS: AsciiLogo[] = [
  {
    logo: `█  █ █▀▀ █▀▀▄ █▀▀▀ █▀▀ █  █ █▀▀█ █▀▀▀
█▀▀█ █▀▀ █  █ █ ▀█ █▀▀ █▀▀█ █  █ █ ▀█
▀  ▀ ▀▀▀ ▀▀▀  ▀▀▀▀ ▀▀▀ ▀  ▀ ▀▀▀▀ ▀▀▀▀`,
  },
  {
    logo: `█▀▀ █▀▀▄ █▀▀▀ █▀▀ █  █
█▀▀ █  █ █ ▀█ █▀▀ █▀▀█
▀▀▀ ▀▀▀  ▀▀▀▀ ▀▀▀ ▀  ▀`,
  },
  {
    logo: 'HEDGEHOG',
  },
];

function getLogoWidth(logo: string): number {
  return logo.split('\n').reduce((max, line) => Math.max(max, line.length), 0);
}

export function selectAsciiLogo(availableWidth: number): AsciiLogo {
  return ASCII_LOGOS.find((entry) => getLogoWidth(entry.logo) <= availableWidth)
    ?? ASCII_LOGOS[ASCII_LOGOS.length - 1];
}
