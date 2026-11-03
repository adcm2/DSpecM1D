(() => {
  const announce = (button, message) => {
    const status = button.nextElementSibling;
    if (status) status.textContent = message;
  };
  document.querySelectorAll('.code-block').forEach((block) => {
    const code = block.querySelector('pre code');
    if (!code) return;
    const button = document.createElement('button');
    button.className = 'copy-code'; button.type = 'button';
    button.setAttribute('aria-label', 'Copy code block');
    button.title = 'Copy code';
    button.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true"><rect x="8" y="8" width="12" height="12" rx="1.5"></rect><path d="M16 8V5.5A1.5 1.5 0 0 0 14.5 4h-9A1.5 1.5 0 0 0 4 5.5v9A1.5 1.5 0 0 0 5.5 16H8"></path></svg>';
    const status = document.createElement('p');
    status.className = 'copy-status'; status.setAttribute('aria-live', 'polite');
    block.insertBefore(button, block.firstChild); block.insertBefore(status, button.nextSibling);
    button.addEventListener('click', async () => {
      const text = code.innerText;
      try {
        if (navigator.clipboard && window.isSecureContext) await navigator.clipboard.writeText(text);
        else {
          const selection = window.getSelection(); const range = document.createRange();
          range.selectNodeContents(code); selection.removeAllRanges(); selection.addRange(range);
          const copied = document.execCommand('copy'); selection.removeAllRanges();
          if (!copied) throw new Error('Copy command unavailable');
        }
        announce(button, 'Copied to clipboard.');
      } catch { announce(button, 'Copy failed. Select the code and copy it manually.'); }
    });
  });
})();
