(() => {
  const currentPath = window.location.pathname;
  document.querySelectorAll("nav a[href]").forEach((link) => {
    const targetPath = new URL(link.href, window.location.href).pathname;
    if (targetPath === currentPath || (currentPath.endsWith("/") && targetPath === `${currentPath}index.html`)) {
      link.setAttribute("aria-current", "page");
    }
  });

  document.querySelectorAll("pre").forEach((pre) => {
    const wrapper = document.createElement("div");
    wrapper.className = "code-block";
    pre.parentNode.insertBefore(wrapper, pre);
    wrapper.appendChild(pre);

    const button = document.createElement("button");
    button.className = "copy-code";
    button.type = "button";
    button.textContent = "Copy code";
    const status = document.createElement("p");
    status.className = "copy-status";
    status.setAttribute("aria-live", "polite");
    status.setAttribute("aria-atomic", "true");
    wrapper.insertBefore(button, pre);
    wrapper.appendChild(status);

    button.addEventListener("click", async () => {
      const text = pre.innerText;
      try {
        if (!navigator.clipboard || !window.isSecureContext) throw new Error("clipboard unavailable");
        await navigator.clipboard.writeText(text);
        status.textContent = "Code copied.";
      } catch (_) {
        status.textContent = "Clipboard access is unavailable. The code is selected; press Ctrl+C or Command+C to copy it.";
        const range = document.createRange();
        range.selectNodeContents(pre);
        const selection = window.getSelection();
        selection.removeAllRanges();
        selection.addRange(range);
        pre.tabIndex = 0;
        pre.focus();
      }
    });
  });
})();
