// This global flag enables manual initialization.
// import CMS, { init } from 'netlify-cms'
const { CMS, initCMS: init } = window;

// Widgets additional
import './app/widgets/index.js';

// Import all shortcodes for editor
import './app/editor/index.js';

// Import collections/config and init
import './app/collections.js';
