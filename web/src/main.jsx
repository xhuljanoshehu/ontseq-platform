import React from "react";
import { createRoot } from "react-dom/client";
import { LiveWorkspace } from "./LiveWorkspace.jsx";
import "./base.css";

createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <LiveWorkspace />
  </React.StrictMode>,
);
