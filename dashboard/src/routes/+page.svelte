<script lang="ts">
  import { invoke } from "@tauri-apps/api/core";
  import { listen, type UnlistenFn } from "@tauri-apps/api/event";
  import { open } from "@tauri-apps/plugin-dialog";
  import { openPath } from "@tauri-apps/plugin-opener";
  import { getCurrentWindow } from "@tauri-apps/api/window";
  import { onMount, onDestroy } from "svelte";
  import { SvelteMap, SvelteSet } from "svelte/reactivity";

  // ── Types ────────────────────────────────────────────────────────────────
  interface Job {
    id: number;
    sample_name: string;
    sample_path: string;
    output_dir: string;
    trees: string;
    status: "running" | "done" | "error";
    start_ts: number;
    end_ts: number | null;
    haplogroup: string | null;
    qc_score: number | null;
    qc1: number | null;
    qc2: number | null;
    qc3: number | null;
    total_reads: number | null;
    valid_markers: number | null;
    reference_genome: string;
    threads: number;
    reads_threshold: number;
    quality_thresh: number;
    base_majority: number;
    prediction_quality: number;
    draw_haplogroups: boolean;
    ancient_dna: boolean;
    private_mutations: boolean;
    collapsed_draw_mode: boolean;
  }

  interface HgPrediction {
    tree: string;
    sample_name: string;
    haplogroup: string;
    hg_marker: string;
    total_reads: number;
    valid_markers: number;
    qc_score: number;
    qc1: number;
    qc2: number;
    qc3: number;
  }

  interface PathNode {
    haplogroup: string;
    marker: string | null;
    state: string;
  }

  interface MarkerStats {
    derived: number;
    ancestral: number;
    not_covered: number;
  }

  interface MixtureContributor {
    rank: number;
    haplogroup: string;
    ratio: number;
    markers: number;
    qc: number;
  }

  interface MixtureSample {
    tree: string;
    sample_name: string;
    common_ancestor: string;
    contributors: MixtureContributor[];
  }

  interface JobResults {
    predictions: HgPrediction[];
    marker_stats: MarkerStats | null;
    mixtures: MixtureSample[];
  }

  // ── Theme ─────────────────────────────────────────────────────────────────
  let darkMode = $state(
    typeof localStorage !== "undefined"
      ? localStorage.getItem("theme") !== "light"
      : true
  );

  $effect(() => {
    document.documentElement.classList.toggle("dark", darkMode);
    localStorage.setItem("theme", darkMode ? "dark" : "light");
  });

  // ── App state ────────────────────────────────────────────────────────────
  type View = "submit" | "result";
  let view = $state<View>("submit");
  let jobs = $state<Job[]>([]);
  let selectedJob = $state<Job | null>(null);
  let jobResults = $state<JobResults | null>(null);
  let pathData = $state<Record<string, PathNode[]> | null>(null);
  let activeTreeTab = $state(0);
  let selectedSample = $state<string | null>(null);
  let resultsError = $state<string | null>(null);

  let ctxMenu = $state<{ x: number; y: number; job: Job } | null>(null);
  let selectedJobIds = $state(new Set<number>());
  let editingJobId = $state<number | null>(null);
  let editingName = $state("");

  function focusAndSelect(node: HTMLInputElement) {
    node.focus();
    node.select();
  }

  function startRename(job: Job) {
    editingJobId = job.id;
    editingName = job.sample_name;
    closeCtxMenu();
  }

  async function saveRename() {
    if (editingJobId === null) return;
    const trimmed = editingName.trim();
    if (trimmed) {
      await invoke("rename_job", { jobId: editingJobId, name: trimmed });
      if (selectedJob?.id === editingJobId) {
        selectedJob = { ...selectedJob, sample_name: trimmed };
      }
      await refreshJobs();
    }
    editingJobId = null;
  }

  function onJobContextMenu(e: MouseEvent, job: Job) {
    e.preventDefault();
    ctxMenu = { x: e.clientX, y: e.clientY, job };
  }

  function closeCtxMenu() { ctxMenu = null; }

  function onJobClick(e: MouseEvent, job: Job) {
    if (e.ctrlKey || e.metaKey) {
      const next = new Set(selectedJobIds);
      if (next.has(job.id)) next.delete(job.id); else next.add(job.id);
      selectedJobIds = next;
    } else {
      selectedJobIds = new Set();
      selectJob(job);
    }
  }

  async function deleteJob(job: Job) {
    await invoke("delete_job", { jobId: job.id });
    if (selectedJob?.id === job.id) {
      selectedJob = null;
      jobResults = null;
      view = "submit";
    }
    jobLogs.delete(job.id);
    runningJobIds.delete(job.id);
    selectedJobIds.delete(job.id);
    await refreshJobs();
    closeCtxMenu();
  }

  async function deleteSelected() {
    if (selectedJobIds.size === 0) return;
    for (const id of selectedJobIds) {
      await invoke("delete_job", { jobId: id });
      if (selectedJob?.id === id) { selectedJob = null; jobResults = null; view = "submit"; }
      jobLogs.delete(id);
      runningJobIds.delete(id);
    }
    selectedJobIds = new Set();
    await refreshJobs();
  }

  function reloadJob(job: Job) {
    bamPath = job.sample_path;
    outputDir = job.output_dir;
    referenceGenome = job.reference_genome;
    selectedTrees = job.trees.split(",");
    threads = job.threads;
    readsThreshold = job.reads_threshold;
    qualityThresh = job.quality_thresh;
    baseMajority = job.base_majority;
    predictionQuality = job.prediction_quality;
    drawHaplogroups = job.draw_haplogroups;
    collapsedDrawMode = job.collapsed_draw_mode;
    ancientDna = job.ancient_dna;
    privateMutations = job.private_mutations;
    selectedJob = null;
    jobResults = null;
    view = "submit";
    closeCtxMenu();
  }

  let sampleNames = $derived([...new Set(jobResults?.predictions.map(p => p.sample_name) ?? [])]);
  let samplePredictions = $derived(
    jobResults?.predictions.filter(p => p.sample_name === selectedSample) ?? []
  );

  // Submit form state
  let bamPath = $state("");
  let outputDir = $state("");
  let referenceGenome = $state("hg38");
  let selectedTrees = $state(["yfull"]);
  let threads = $state(4);
  let readsThreshold = $state(10);
  let qualityThresh = $state(20);
  let baseMajority = $state(90);
  let predictionQuality = $state(0.95);
  let drawHaplogroups = $state(false);
  let collapsedDrawMode = $state(false);
  let ancientDna = $state(false);
  let privateMutations = $state(false);
  let mixtureMode = $state(false);

  interface ProgressBar { stage: string; current: number; total: number; }
  interface ProgressBarPayload { job_id: number; stage: string; current: number; total: number; }

  // Multi-run state
  let runningJobIds = new SvelteSet<number>();
  let jobLogs = new SvelteMap<number, string[]>();
  let jobProgress = new SvelteMap<number, ProgressBar[]>();
  let logJobId = $state<number | null>(null);
  let logEl: HTMLElement | null = null;
  let unlisteners: UnlistenFn[] = [];

  let displayLog = $derived(logJobId !== null ? (jobLogs.get(logJobId) ?? []) : []);

  const TREES = [
    { id: "yfull",     label: "YFull v14" },
    { id: "yfull_v10", label: "YFull v10" },
    { id: "ftdna",     label: "FTDNA" },
    { id: "isogg",     label: "ISOGG" },
  ];

  const REF_GENOMES = [
    { id: "hg38", label: "hg38" },
    { id: "hg19", label: "hg19" },
    { id: "t2t",  label: "T2T" },
  ];

  // ── Helpers ──────────────────────────────────────────────────────────────
  $effect(() => {
    displayLog.length;
    if (logEl) logEl.scrollTop = logEl.scrollHeight;
  });

  function fmt_ts(ts: number): string {
    return new Date(ts * 1000).toLocaleString();
  }

  function qc_color(v: number): string {
    if (v >= 0.95) return "#6bcb77";
    if (v >= 0.8)  return "#f9c74f";
    return "#e94560";
  }

  function pct(n: number, total: number): string {
    return total > 0 ? ((n / total) * 100).toFixed(1) + "%" : "0%";
  }

  // ── Job list ─────────────────────────────────────────────────────────────
  async function refreshJobs() {
    jobs = await invoke<Job[]>("list_jobs");
  }

  async function selectJob(job: Job) {
    selectedJob = job;
    jobResults = null;
    pathData = null;
    resultsError = null;
    activeTreeTab = 0;
    selectedSample = null;
    logJobId = job.id;
    view = "result";
    if (job.status === "done") {
      try {
        jobResults = await invoke<JobResults>("get_job_results", { outputDir: job.output_dir, trees: job.trees });
        selectedSample = jobResults?.predictions[0]?.sample_name ?? null;
      } catch (e) {
        resultsError = String(e);
      }
    }
  }

  // Fetch per-sample path data whenever the selected sample changes
  $effect(() => {
    if (selectedJob?.status === "done" && selectedSample) {
      const trees = selectedJob.trees.split(",");
      invoke<Record<string, PathNode[]>>("get_haplogroup_path", {
        outputDir: selectedJob.output_dir,
        sampleName: selectedSample,
        trees,
      }).then(data => { pathData = data; });
    } else {
      pathData = null;
    }
  });

  function treeHtmlPaths(job: Job): string[] {
    const trees = job.trees.split(",");
    if (trees.length === 1) {
      return [`${job.output_dir}/hg_tree_image.html`];
    }
    return trees.map(t => `${job.output_dir}/hg_tree_image_${t}.html`);
  }

  async function openTreeHtml(path: string) {
    try { await openPath(path); }
    catch (e) { alert(`Could not open ${path}\n${e}`); }
  }

  function newRun() {
    selectedJob = null;
    jobResults = null;
    view = "submit";
  }

  // ── Run ──────────────────────────────────────────────────────────────────
  const INPUT_FILTERS = [
    { name: "All supported",    extensions: ["bam", "cram", "vcf", "vcf.gz", "fastq", "fastq.gz", "fq", "fq.gz", "ped", "bed"] },
    { name: "BAM / CRAM",       extensions: ["bam", "cram"] },
    { name: "VCF",              extensions: ["vcf", "vcf.gz"] },
    { name: "FASTQ",            extensions: ["fastq", "fastq.gz", "fq", "fq.gz"] },
    { name: "PLINK",            extensions: ["ped", "bed"] },
    { name: "All files",        extensions: ["*"] },
  ];

  async function browseInputFile() {
    try { await getCurrentWindow().maximize(); } catch (_) {}
    const result = await open({ multiple: false, directory: false, filters: INPUT_FILTERS });
    if (result) bamPath = result as string;
  }

  async function browseInputDir() {
    try { await getCurrentWindow().maximize(); } catch (_) {}
    const result = await open({ multiple: false, directory: true });
    if (result) bamPath = result as string;
  }

  async function browseOutputDir() {
    try { await getCurrentWindow().maximize(); } catch (_) {}
    const result = await open({ multiple: false, directory: true });
    if (result) outputDir = result as string;
  }

  function toggleTree(id: string) {
    selectedTrees = selectedTrees.includes(id)
      ? selectedTrees.filter((t) => t !== id)
      : [...selectedTrees, id];
  }

  let isLaunching = $state(false);

  async function run() {
    if (!bamPath.trim() || !outputDir.trim() || selectedTrees.length === 0) return;
    if (isLaunching) return;
    isLaunching = true;
    try {
      const jobId = await invoke<number>("run_yleaf", {
        bamPath: bamPath.trim(),
        outputDir: outputDir.trim(),
        referenceGenome,
        tree: selectedTrees,
        threads,
        readsThreshold,
        qualityThresh,
        baseMajority,
        predictionQuality,
        drawHaplogroups,
        collapsedDrawMode,
        ancientDna,
        privateMutations,
        mixtureMode,
      });
      runningJobIds.add(jobId);
      jobLogs.set(jobId, []);
      logJobId = jobId;
      await refreshJobs();
    } catch (e) {
      jobLogs.set(-1, [`Failed to start: ${String(e)}`]);
      logJobId = -1;
    } finally {
      isLaunching = false;
    }
  }

  async function stopRun(jobId: number) {
    try {
      await invoke("kill_job", { jobId });
      const lines = jobLogs.get(jobId) ?? [];
      jobLogs.set(jobId, [...lines, "", "Run stopped by user."]);
    } catch (e) {
      const lines = jobLogs.get(jobId) ?? [];
      jobLogs.set(jobId, [...lines, "", `Stop failed: ${String(e)}`]);
    }
    runningJobIds.delete(jobId);
    await refreshJobs();
  }

  let canRun = $derived(
    bamPath.trim() !== "" &&
    outputDir.trim() !== "" &&
    selectedTrees.length > 0 &&
    !isLaunching
  );

  let _onKeyDown: ((e: KeyboardEvent) => void) | null = null;

  onMount(async () => {
    await refreshJobs();
    _onKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Delete" && selectedJobIds.size > 0) deleteSelected();
    };
    window.addEventListener("keydown", _onKeyDown);
    unlisteners = await Promise.all([
      listen<{ job_id: number; line: string }>("yleaf-progress", (e) => {
        const { job_id, line } = e.payload;
        const lines = jobLogs.get(job_id) ?? [];
        lines.push(line);
        jobLogs.set(job_id, lines);
      }),
      listen<ProgressBarPayload>("yleaf-progress-bar", (e) => {
        const { job_id, stage, current, total } = e.payload;
        const bars = [...(jobProgress.get(job_id) ?? [])];
        const idx = bars.findIndex(b => b.stage === stage);
        if (idx >= 0) bars[idx] = { stage, current, total };
        else bars.push({ stage, current, total });
        jobProgress.set(job_id, bars);
      }),
      listen<{ job_id: number; output_dir: string }>("yleaf-done", async (e) => {
        const { job_id, output_dir } = e.payload;
        jobLogs.set(job_id, [...(jobLogs.get(job_id) ?? []), "", `Done. Output: ${output_dir}`]);
        runningJobIds.delete(job_id);
        await refreshJobs();
      }),
      listen<[number, string]>("yleaf-error", async (e) => {
        const [job_id, msg] = e.payload;
        jobLogs.set(job_id, [...(jobLogs.get(job_id) ?? []), "", `Error: ${msg}`]);
        runningJobIds.delete(job_id);
        await refreshJobs();
      }),
      listen<number>("job-updated", async () => {
        await refreshJobs();
      }),
    ]);
  });
  onDestroy(() => {
    for (const u of unlisteners) u();
    if (_onKeyDown) window.removeEventListener("keydown", _onKeyDown);
  });
</script>

<!-- ══════════════════════════════════════════════════════════════════════ -->
<div class="flex h-screen overflow-hidden bg-gray-100 dark:bg-navy text-slate-800 dark:text-pale">

  <!-- ── Sidebar ── -->
  <nav class="w-52 shrink-0 flex flex-col overflow-y-auto
              bg-white dark:bg-ink
              border-r border-slate-200 dark:border-well">

    <!-- Header -->
    <div class="flex items-center gap-1.5 px-3 py-2.5 shrink-0
                border-b border-slate-200 dark:border-well">
      {#if selectedJobIds.size > 0}
        <span class="flex-1 text-[0.62rem] font-medium uppercase tracking-widest text-red-500 dark:text-brand">
          {selectedJobIds.size} selected
        </span>
        <button
          onclick={deleteSelected}
          class="text-[0.65rem] px-2 py-0.5 rounded border-none cursor-pointer
                 bg-red-500 hover:opacity-85 transition-opacity text-white font-medium">
          Delete
        </button>
      {:else}
        <span class="flex-1 text-[0.62rem] font-medium uppercase tracking-widest
                     text-slate-400 dark:text-muted">
          History
        </span>
      {/if}
      <!-- Theme toggle -->
      <button
        onclick={() => (darkMode = !darkMode)}
        title={darkMode ? "Switch to light mode" : "Switch to dark mode"}
        class="p-1 rounded transition-colors
               text-slate-400 hover:text-slate-600
               dark:text-muted dark:hover:text-pale">
        {#if darkMode}
          <!-- Sun: switch to light -->
          <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round">
            <circle cx="12" cy="12" r="5"/>
            <line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/>
            <line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/>
            <line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/>
            <line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/>
          </svg>
        {:else}
          <!-- Moon: switch to dark -->
          <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round">
            <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/>
          </svg>
        {/if}
      </button>
      <button
        onclick={newRun}
        class="text-[0.65rem] px-2 py-0.5 rounded border-none cursor-pointer
               bg-brand hover:opacity-85 transition-opacity text-white font-medium">
        + New
      </button>
    </div>

    <!-- Job list -->
    {#if jobs.length === 0}
      <p class="m-0 px-3 py-3 text-[0.73rem] text-slate-400 dark:text-ghost">No runs yet.</p>
    {/if}
    {#each jobs as job}
      <!-- svelte-ignore a11y_no_static_element_interactions -->
      <div
        class="border-b px-3 py-2 cursor-pointer flex flex-col gap-0.5 transition-colors
               border-slate-100 dark:border-[#111]
               {selectedJobIds.has(job.id)
                 ? 'bg-red-50 dark:bg-[#2a1010]'
                 : selectedJob?.id === job.id
                   ? 'bg-blue-50 dark:bg-raised'
                   : 'hover:bg-slate-50 dark:hover:bg-panel'}"
        onclick={(e) => onJobClick(e, job)}
        onkeydown={(e) => e.key === "Enter" && selectJob(job)}
        oncontextmenu={(e) => onJobContextMenu(e, job)}
        role="button"
        tabindex="0"
      >
        {#if editingJobId === job.id}
          <!-- svelte-ignore a11y_autofocus -->
          <input
            type="text"
            bind:value={editingName}
            use:focusAndSelect
            onkeydown={(e) => {
              if (e.key === "Enter") { e.preventDefault(); saveRename(); }
              if (e.key === "Escape") { editingJobId = null; }
              e.stopPropagation();
            }}
            onclick={(e) => e.stopPropagation()}
            onblur={saveRename}
            class="text-[0.78rem] font-semibold w-full bg-transparent outline-none
                   border-b border-blue-400 dark:border-azure
                   text-slate-800 dark:text-pale"
          />
        {:else}
          <span class="text-[0.78rem] font-semibold truncate
                       text-slate-800 dark:text-pale">
            {job.sample_name}
          </span>
        {/if}
        <span class="text-[0.7rem] truncate text-sky-600 dark:text-teal">
          {job.haplogroup ?? "—"}
        </span>
        <div class="flex items-center gap-1 flex-wrap">
          <span class="text-[0.58rem] px-1.5 py-0.5 rounded-full font-bold uppercase
            {job.status === 'running' ? 'bg-blue-100 text-blue-700 dark:bg-[#1a4a8a] dark:text-teal'
           : job.status === 'done'    ? 'bg-green-100 text-green-700 dark:bg-[#1a3a2a] dark:text-lime'
           :                            'bg-red-100 text-red-600 dark:bg-[#3a1a1a] dark:text-brand'}">
            {job.status}
          </span>
          <span class="text-[0.62rem] text-slate-400 dark:text-ghost">{fmt_ts(job.start_ts)}</span>
          {#if runningJobIds.has(job.id)}
            <button
              class="ml-auto border-none rounded text-[0.58rem] px-1.5 py-0.5 cursor-pointer leading-none transition-colors
                     bg-red-100 text-red-600 hover:bg-red-200
                     dark:bg-[#5a2020] dark:text-brand dark:hover:bg-[#7a2020]"
              onclick={(e) => { e.stopPropagation(); stopRun(job.id); }}
              title="Stop run">&#9632;</button>
          {/if}
        </div>
      </div>
    {/each}
  </nav>

  <!-- ── Center panel ── -->
  <main class="flex-1 overflow-y-auto border-r border-slate-200 dark:border-well">

    {#if view === "submit"}
      <div class="p-6 flex flex-col gap-5 w-full">

        <!-- Brand header -->
        <div class="flex items-center gap-4 pb-5 border-b border-slate-200 dark:border-well">
          <img src="/yleaf_logo.png" alt="Yleaf" class="h-40 w-auto" />
          <div class="flex-1">
            <h1 class="m-0 text-2xl font-bold text-brand tracking-[0.12em]">Yleaf 4.0</h1>
            <p class="m-0 mt-1 text-xs text-slate-500 dark:text-muted">
              Y-chromosome haplogroup inference
            </p>
          </div>
          <img src="/Erasmus-MC-logo.png" alt="Erasmus MC" class="h-16 w-auto opacity-80 ml-auto dark:opacity-90" />
        </div>

        <!-- Input -->
        <div class="flex flex-col gap-2">
          <span class="text-[0.62rem] font-medium uppercase tracking-widest
                       text-slate-500 dark:text-muted">
            Input &mdash; BAM / CRAM / VCF / FASTQ / PLINK / directory
          </span>
          <div class="rounded-lg px-3 py-2 text-[0.76rem] break-all min-h-[2.2rem] leading-relaxed border
                      bg-slate-50 border-slate-200 dark:bg-well dark:border-rim
                      {bamPath ? 'text-slate-800 dark:text-pale' : 'text-slate-400 dark:text-ghost'}">
            {bamPath || "No file selected"}
          </div>
          <div class="flex gap-2">
            <button onclick={browseInputFile}
              class="flex-1 rounded-lg border py-1.5 text-[0.75rem] font-medium cursor-pointer transition-colors
                     bg-white border-slate-300 text-slate-700 hover:bg-slate-50
                     dark:bg-rim dark:border-[#2a6aaa] dark:text-pale dark:hover:bg-[#2a5a9a]">
              Browse file
            </button>
            <button onclick={browseInputDir}
              class="flex-1 rounded-lg border py-1.5 text-[0.75rem] font-medium cursor-pointer transition-colors
                     bg-white border-slate-300 text-slate-700 hover:bg-slate-50
                     dark:bg-rim dark:border-[#2a6aaa] dark:text-pale dark:hover:bg-[#2a5a9a]">
              Browse folder
            </button>
          </div>
        </div>

        <!-- Output -->
        <div class="flex flex-col gap-2">
          <span class="text-[0.62rem] font-medium uppercase tracking-widest
                       text-slate-500 dark:text-muted">
            Output directory
          </span>
          <input
            type="text"
            placeholder="Type a path or browse…"
            bind:value={outputDir}
            class="w-full rounded-lg border px-3 py-2 text-[0.82rem] outline-none transition-colors
                   bg-white border-slate-300 text-slate-800 placeholder:text-slate-400 focus:border-blue-400
                   dark:bg-well dark:border-rim dark:text-pale dark:placeholder:text-ghost dark:focus:border-azure"
          />
          <button onclick={browseOutputDir}
            class="self-start rounded-lg border px-3 py-1.5 text-[0.75rem] font-medium cursor-pointer transition-colors
                   bg-white border-slate-300 text-slate-700 hover:bg-slate-50
                   dark:bg-rim dark:border-[#2a6aaa] dark:text-pale dark:hover:bg-[#2a5a9a]">
            Browse folder
          </button>
        </div>

        <!-- Reference genome -->
        <div class="flex flex-col gap-2">
          <span class="text-[0.62rem] font-medium uppercase tracking-widest
                       text-slate-500 dark:text-muted">
            Reference genome
          </span>
          <div class="flex gap-5">
            {#each REF_GENOMES as rg}
              <label class="flex items-center gap-2 text-[0.82rem] cursor-pointer
                            text-slate-700 dark:text-pale">
                <input type="radio" name="ref" value={rg.id} bind:group={referenceGenome}
                  class="accent-brand" />
                {rg.label}
              </label>
            {/each}
          </div>
        </div>

        <!-- Trees -->
        <div class="flex flex-col gap-2">
          <span class="text-[0.62rem] font-medium uppercase tracking-widest
                       text-slate-500 dark:text-muted">
            Trees
          </span>
          <div class="flex gap-5 flex-wrap">
            {#each TREES as tree}
              <label class="flex items-center gap-2 text-[0.82rem] cursor-pointer
                            text-slate-700 dark:text-pale">
                <input
                  type="checkbox"
                  checked={selectedTrees.includes(tree.id)}
                  onchange={() => toggleTree(tree.id)}
                  class="accent-brand"
                />
                {tree.label}
              </label>
            {/each}
          </div>
        </div>

        <!-- Advanced options -->
        <div class="flex flex-col gap-3">
          <span class="text-[0.62rem] font-medium uppercase tracking-widest
                       text-slate-500 dark:text-muted">
            Advanced options
          </span>
          <div class="grid gap-x-6 gap-y-3 pl-1"
               style="grid-template-columns: repeat(auto-fill, minmax(160px, 1fr))">
            <label class="flex flex-col gap-1.5">
              <span class="text-[0.62rem] font-medium uppercase tracking-widest
                           text-slate-500 dark:text-muted">Threads (-t)</span>
              <input type="number" min="1" max="64" bind:value={threads}
                class="rounded-lg border px-2.5 py-1.5 text-[0.82rem] outline-none w-full transition-colors
                       bg-white border-slate-300 text-slate-800 focus:border-blue-400
                       dark:bg-well dark:border-rim dark:text-pale dark:focus:border-azure" />
            </label>
            <label class="flex flex-col gap-1.5">
              <span class="text-[0.62rem] font-medium uppercase tracking-widest
                           text-slate-500 dark:text-muted">Min reads (-r)</span>
              <input type="number" min="1" bind:value={readsThreshold}
                class="rounded-lg border px-2.5 py-1.5 text-[0.82rem] outline-none w-full transition-colors
                       bg-white border-slate-300 text-slate-800 focus:border-blue-400
                       dark:bg-well dark:border-rim dark:text-pale dark:focus:border-azure" />
            </label>
            <label class="flex flex-col gap-1.5">
              <span class="text-[0.62rem] font-medium uppercase tracking-widest
                           text-slate-500 dark:text-muted">Read quality (-q)</span>
              <input type="number" min="10" max="40" bind:value={qualityThresh}
                class="rounded-lg border px-2.5 py-1.5 text-[0.82rem] outline-none w-full transition-colors
                       bg-white border-slate-300 text-slate-800 focus:border-blue-400
                       dark:bg-well dark:border-rim dark:text-pale dark:focus:border-azure" />
            </label>
            <label class="flex flex-col gap-1.5">
              <span class="text-[0.62rem] font-medium uppercase tracking-widest
                           text-slate-500 dark:text-muted">Base majority % (-b)</span>
              <input type="number" min="50" max="99" bind:value={baseMajority}
                class="rounded-lg border px-2.5 py-1.5 text-[0.82rem] outline-none w-full transition-colors
                       bg-white border-slate-300 text-slate-800 focus:border-blue-400
                       dark:bg-well dark:border-rim dark:text-pale dark:focus:border-azure" />
            </label>
            <label class="flex flex-col gap-1.5">
              <span class="text-[0.62rem] font-medium uppercase tracking-widest
                           text-slate-500 dark:text-muted">Pred. quality (-pq)</span>
              <input type="number" min="0" max="1" step="0.01" bind:value={predictionQuality}
                class="rounded-lg border px-2.5 py-1.5 text-[0.82rem] outline-none w-full transition-colors
                       bg-white border-slate-300 text-slate-800 focus:border-blue-400
                       dark:bg-well dark:border-rim dark:text-pale dark:focus:border-azure" />
            </label>
            <div class="col-span-full flex flex-col gap-2 pt-1">
              <label class="flex items-center gap-2 text-[0.82rem] cursor-pointer
                            text-slate-700 dark:text-pale">
                <input type="checkbox" bind:checked={drawHaplogroups} class="accent-brand" />
                Draw haplogroup tree (-dh)
              </label>
              <label class="flex items-center gap-2 text-[0.82rem] cursor-pointer
                            {drawHaplogroups ? 'text-slate-700 dark:text-pale' : 'text-slate-400 dark:text-ghost'}">
                <input type="checkbox" bind:checked={collapsedDrawMode} disabled={!drawHaplogroups} class="accent-brand" />
                Collapse tree image (-hc)
              </label>
              <label class="flex items-center gap-2 text-[0.82rem] cursor-pointer
                            text-slate-700 dark:text-pale">
                <input type="checkbox" bind:checked={ancientDna} class="accent-brand" />
                Ancient DNA mode (-aDNA)
              </label>
              <label class="flex items-center gap-2 text-[0.82rem] cursor-pointer
                            text-slate-700 dark:text-pale">
                <input type="checkbox" bind:checked={privateMutations} class="accent-brand" />
                Find private mutations (-p)
              </label>
              <label class="flex items-center gap-2 text-[0.82rem] cursor-pointer
                            text-slate-700 dark:text-pale">
                <input type="checkbox" bind:checked={mixtureMode} class="accent-brand" />
                Mixture mode (-mix)
              </label>
            </div>
          </div>
        </div>

        <!-- Run button -->
        <button
          onclick={run}
          disabled={!canRun}
          class="w-full py-3 rounded-xl border-none font-semibold text-sm cursor-pointer transition-opacity
                 bg-brand text-white
                 hover:opacity-90 disabled:opacity-40 disabled:cursor-not-allowed">
          Run Yleaf
        </button>

      </div>

    {:else if view === "result" && selectedJob}
      <div class="p-6 flex flex-col gap-5">

        <!-- Result header -->
        <div class="flex items-start justify-between gap-3">
          <div class="min-w-0">
            <h2 class="m-0 text-lg font-semibold text-slate-800 dark:text-pale">
              {selectedJob.sample_name}
            </h2>
            <span class="text-[0.68rem] break-all text-slate-400 dark:text-ghost">
              {selectedJob.sample_path}
            </span>
          </div>
          <span class="shrink-0 text-[0.58rem] px-2 py-0.5 rounded-full font-bold uppercase
            {selectedJob.status === 'running' ? 'bg-blue-100 text-blue-700 dark:bg-[#1a4a8a] dark:text-teal'
           : selectedJob.status === 'done'    ? 'bg-green-100 text-green-700 dark:bg-[#1a3a2a] dark:text-lime'
           :                                    'bg-red-100 text-red-600 dark:bg-[#3a1a1a] dark:text-brand'}">
            {selectedJob.status}
          </span>
        </div>

        {#if selectedJob.status === "running"}
          <p class="m-0 text-[0.82rem] text-slate-500 dark:text-muted">
            Run in progress — results will appear when complete.
          </p>

        {:else if selectedJob.status === "error"}
          <p class="m-0 text-[0.82rem] text-red-500 dark:text-brand">
            Run failed. Check log for details.
          </p>

        {:else if resultsError}
          <!-- Output files missing — show cached DB values -->
          <p class="m-0 text-[0.76rem] text-amber-600 dark:text-yellow-400">
            Output files unavailable (may have been moved or deleted). Showing cached results.
          </p>
          {#if selectedJob.haplogroup}
            <div class="rounded-xl p-4 border bg-white border-slate-200 dark:bg-panel dark:border-well">
              <div class="text-2xl font-bold mb-1.5 text-sky-600 dark:text-teal">
                {selectedJob.haplogroup}
              </div>
              <div class="text-[0.73rem] text-slate-400 dark:text-ghost">
                {(selectedJob.total_reads ?? 0).toLocaleString()} mapped reads &middot; {(selectedJob.valid_markers ?? 0).toLocaleString()} markers
              </div>
            </div>
            {#if selectedJob.qc_score !== null}
              <div class="grid grid-cols-2 gap-2">
                {#each [
                  { label: "Total QC",        v: selectedJob.qc_score },
                  { label: "QC-1 (backbone)", v: selectedJob.qc1 ?? 0 },
                  { label: "QC-2 (coverage)", v: selectedJob.qc2 ?? 0 },
                  { label: "QC-3 (data)",     v: selectedJob.qc3 ?? 0 },
                ] as q}
                  <div class="rounded-xl p-3 border bg-slate-50 border-slate-200 dark:bg-void dark:border-transparent">
                    <span class="block text-[0.6rem] font-medium uppercase tracking-widest mb-1.5 text-slate-400 dark:text-muted">{q.label}</span>
                    <span class="block text-base font-bold mb-2" style="color:{qc_color(q.v)}">{q.v.toFixed(3)}</span>
                    <div class="rounded-full h-1 bg-slate-200 dark:bg-navy">
                      <div class="h-1 rounded-full" style="width:{(q.v*100).toFixed(1)}%;background:{qc_color(q.v)}"></div>
                    </div>
                  </div>
                {/each}
              </div>
            {/if}
          {/if}
          <!-- Run metadata still shown -->
          <div class="flex flex-col gap-1 border-t pt-4 border-slate-200 dark:border-well">
            <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Started: {fmt_ts(selectedJob.start_ts)}</span>
            {#if selectedJob.end_ts}
              <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Finished: {fmt_ts(selectedJob.end_ts)}</span>
              <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Duration: {Math.round((selectedJob.end_ts - selectedJob.start_ts) / 60)} min</span>
            {/if}
            <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Trees: {selectedJob.trees}</span>
            <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Output: {selectedJob.output_dir}</span>
          </div>

        {:else if jobResults}
          <!-- Sample selector -->
          {#if sampleNames.length > 1}
            <div class="flex items-center gap-2 text-[0.82rem] text-slate-500 dark:text-muted">
              <label for="sample-select" class="shrink-0">Sample</label>
              <select id="sample-select" bind:value={selectedSample}
                onchange={() => { activeTreeTab = 0; }}
                class="flex-1 max-w-xs rounded-lg border px-2.5 py-1.5 text-[0.82rem] cursor-pointer outline-none
                       bg-white border-slate-300 text-slate-800
                       dark:bg-[#0f1a2e] dark:border-[#1e3a6e] dark:text-pale">
                {#each sampleNames as name}
                  <option value={name}>{name}</option>
                {/each}
              </select>
            </div>
            {#if selectedSample}
              <h3 class="m-0 text-base font-semibold text-slate-700 dark:text-pale">
                {selectedSample}
              </h3>
            {/if}
          {/if}

          <!-- Tree tabs -->
          {#if samplePredictions.length > 1}
            <div class="flex gap-1 flex-wrap">
              {#each samplePredictions as pred, i}
                <button
                  onclick={() => (activeTreeTab = i)}
                  class="px-3 py-1.5 rounded-t-lg border text-[0.76rem] cursor-pointer transition-colors
                    {activeTreeTab === i
                      ? 'border-slate-300 bg-white text-slate-800 dark:bg-panel dark:border-rim dark:text-pale'
                      : 'border-slate-200 bg-slate-100 text-slate-500 hover:text-slate-700 dark:bg-well dark:border-rim dark:text-muted dark:hover:text-pale'}">
                  {pred.tree}
                </button>
              {/each}
            </div>
          {/if}

          {#each samplePredictions as pred, i}
            {#if i === activeTreeTab}
              <!-- Haplogroup card -->
              <div class="rounded-xl p-4 border
                          bg-white border-slate-200
                          dark:bg-panel dark:border-well">
                <div class="text-2xl font-bold mb-1.5 text-sky-600 dark:text-teal">
                  {pred.haplogroup}
                </div>
                <div class="text-[0.82rem] mb-1 text-slate-500 dark:text-muted">
                  Defining marker: <strong class="text-slate-700 dark:text-pale">{pred.hg_marker}</strong>
                </div>
                <div class="text-[0.73rem] text-slate-400 dark:text-ghost">
                  {pred.total_reads.toLocaleString()} mapped reads &middot; {pred.valid_markers.toLocaleString()} markers
                </div>
              </div>

              <!-- QC scores -->
              <div class="grid grid-cols-2 gap-2">
                {#each [
                  { label: "Total QC",       v: pred.qc_score },
                  { label: "QC-1 (backbone)", v: pred.qc1 },
                  { label: "QC-2 (coverage)", v: pred.qc2 },
                  { label: "QC-3 (data)",     v: pred.qc3 },
                ] as q}
                  <div class="rounded-xl p-3 border
                              bg-slate-50 border-slate-200
                              dark:bg-void dark:border-transparent">
                    <span class="block text-[0.6rem] font-medium uppercase tracking-widest mb-1.5
                                 text-slate-400 dark:text-muted">
                      {q.label}
                    </span>
                    <span class="block text-base font-bold mb-2" style="color:{qc_color(q.v)}">
                      {q.v.toFixed(3)}
                    </span>
                    <div class="rounded-full h-1 bg-slate-200 dark:bg-navy">
                      <div class="h-1 rounded-full transition-all"
                           style="width:{(q.v * 100).toFixed(1)}%;background:{qc_color(q.v)}">
                      </div>
                    </div>
                  </div>
                {/each}
              </div>
            {/if}
          {/each}

          <!-- Haplogroup path — per selected sample -->
          {#each samplePredictions as pred, i}
            {#if i === activeTreeTab && pathData?.[pred.tree]}
              {@const measured = pathData[pred.tree].filter(n => n.state === "D" || n.state === "A")}
              {#if measured.length > 0}
                <div class="border-t pt-4 border-slate-200 dark:border-well">
                  <h3 class="m-0 mb-3 text-[0.65rem] font-semibold uppercase tracking-widest
                             text-slate-400 dark:text-muted">
                    Haplogroup path
                  </h3>
                  <div class="flex flex-col pl-1">
                    {#each measured as node, ni}
                      <div class="flex items-start relative">
                        <div class="w-2.5 h-2.5 rounded-full mt-0.5 shrink-0"
                             style="background:{node.state === 'D' ? '#6bcb77' : '#e94560'}">
                        </div>
                        <div class="ml-2 flex items-center gap-2 pb-1">
                          <span class="text-[0.76rem] font-semibold
                                       text-slate-700 dark:text-[#c8d6e5]">
                            {node.haplogroup}
                          </span>
                          {#if node.marker}
                            <span class="text-[0.65rem] font-mono px-1 py-px rounded
                              {node.state === 'D'
                                ? 'bg-green-100 text-green-700 dark:bg-[rgba(107,203,119,0.15)] dark:text-lime'
                                : 'bg-red-100 text-red-600 dark:bg-[rgba(233,69,96,0.15)] dark:text-brand'}">
                              {node.marker} [{node.state}]
                            </span>
                          {/if}
                        </div>
                        {#if ni < measured.length - 1}
                          <div class="absolute left-[4px] top-[13px] w-0.5 bg-slate-200 dark:bg-well"
                               style="height:calc(100% - 2px)">
                          </div>
                        {/if}
                      </div>
                    {/each}
                  </div>
                </div>
              {/if}
            {/if}
          {/each}

          <!-- Marker stats -->
          {#if jobResults.marker_stats}
            {@const ms = jobResults.marker_stats}
            {@const total = ms.derived + ms.ancestral + ms.not_covered}
            <div class="border-t pt-4 border-slate-200 dark:border-well">
              <h3 class="m-0 mb-3 text-[0.65rem] font-semibold uppercase tracking-widest
                         text-slate-400 dark:text-muted">
                Marker summary
              </h3>
              <div class="flex flex-col gap-2">
                {#each [
                  { label: "Derived (D)",  n: ms.derived,     color: "#6bcb77" },
                  { label: "Ancestral (A)", n: ms.ancestral,  color: "#4a9eff" },
                  { label: "Not covered",  n: ms.not_covered, color: "#888"    },
                ] as m}
                  <div class="grid items-center gap-2"
                       style="grid-template-columns: 110px 1fr 130px">
                    <span class="text-[0.73rem] text-slate-500 dark:text-muted">{m.label}</span>
                    <div class="rounded-full h-1.5 bg-slate-200 dark:bg-navy">
                      <div class="h-1.5 rounded-full transition-all"
                           style="width:{pct(m.n,total)};background:{m.color}">
                      </div>
                    </div>
                    <span class="text-[0.7rem] text-right text-slate-400 dark:text-ghost">
                      {m.n.toLocaleString()} ({pct(m.n, total)})
                    </span>
                  </div>
                {/each}
              </div>
            </div>
          {/if}

          <!-- Mixture results -->
          {#if jobResults.mixtures && jobResults.mixtures.filter(m => m.sample_name === selectedSample && m.tree === samplePredictions[activeTreeTab]?.tree).length > 0}
            <div class="border-t pt-4 border-slate-200 dark:border-well">
              <h3 class="m-0 mb-3 text-[0.65rem] font-semibold uppercase tracking-widest
                         text-slate-400 dark:text-muted">
                Mixture analysis
              </h3>
              {#each jobResults.mixtures.filter(m => m.sample_name === selectedSample && m.tree === samplePredictions[activeTreeTab]?.tree) as mix}
                <div class="mb-3 rounded border border-amber-300 dark:border-amber-700
                            bg-amber-50 dark:bg-amber-950/30 p-3">
                  <div class="flex items-center gap-2 mb-2">
                    <span class="rounded bg-amber-400 dark:bg-amber-600 px-1.5 py-0.5
                                 text-[0.6rem] font-bold uppercase tracking-wider text-white">
                      MIX
                    </span>
                    <span class="text-[0.72rem] text-slate-500 dark:text-muted">
                      Tree: {mix.tree} · LCA: {mix.common_ancestor}
                      · {mix.contributors.length} contributors
                    </span>
                  </div>
                  <table class="w-full text-[0.72rem] border-collapse">
                    <thead>
                      <tr class="text-slate-400 dark:text-muted text-left">
                        <th class="pb-1 pr-3 font-medium">#</th>
                        <th class="pb-1 pr-3 font-medium">Haplogroup</th>
                        <th class="pb-1 pr-3 font-medium">Ratio</th>
                        <th class="pb-1 pr-3 font-medium">Markers</th>
                        <th class="pb-1 font-medium">QC</th>
                      </tr>
                    </thead>
                    <tbody>
                      {#each mix.contributors as c}
                        <tr class="text-slate-700 dark:text-pale border-t
                                   border-slate-100 dark:border-well">
                          <td class="py-0.5 pr-3">{c.rank}</td>
                          <td class="py-0.5 pr-3 font-mono">{c.haplogroup}</td>
                          <td class="py-0.5 pr-3">{(c.ratio * 100).toFixed(1)}%</td>
                          <td class="py-0.5 pr-3">{c.markers}</td>
                          <td class="py-0.5">{c.qc.toFixed(3)}</td>
                        </tr>
                      {/each}
                    </tbody>
                  </table>
                </div>
              {/each}
            </div>
          {/if}

          <!-- Run metadata -->
          <div class="flex flex-col gap-1 border-t pt-4 border-slate-200 dark:border-well">
            <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Started: {fmt_ts(selectedJob.start_ts)}</span>
            {#if selectedJob.end_ts}
              <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Finished: {fmt_ts(selectedJob.end_ts)}</span>
              <span class="text-[0.7rem] text-slate-400 dark:text-ghost">
                Duration: {Math.round((selectedJob.end_ts - selectedJob.start_ts) / 60)} min
              </span>
            {/if}
            <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Trees: {selectedJob.trees}</span>
            <span class="text-[0.7rem] text-slate-400 dark:text-ghost">Output: {selectedJob.output_dir}</span>
            {#if selectedJob.draw_haplogroups}
              <div class="flex flex-wrap gap-2 pt-1">
                {#each treeHtmlPaths(selectedJob) as htmlPath}
                  <button
                    onclick={() => openTreeHtml(htmlPath)}
                    class="rounded-lg border px-3 py-1.5 text-[0.75rem] font-medium cursor-pointer transition-colors
                           bg-white border-slate-300 text-slate-700 hover:bg-slate-50
                           dark:bg-rim dark:border-[#2a6aaa] dark:text-pale dark:hover:bg-[#2a5a9a]">
                    Open tree — {htmlPath.split('/').pop()}
                  </button>
                {/each}
              </div>
            {/if}
          </div>
        {/if}
      </div>
    {/if}
  </main>

  <!-- ── Log panel ── -->
  <div class="w-96 shrink-0 flex flex-col p-3 overflow-hidden
              bg-slate-50 dark:bg-void">
    <div class="flex items-center justify-between pb-2 mb-2 shrink-0
                border-b border-slate-200 dark:border-well">
      <span class="text-[0.62rem] font-medium uppercase tracking-widest
                   text-slate-400 dark:text-muted">
        Log{logJobId !== null ? ` · job #${logJobId}` : ""}
      </span>
      {#if displayLog.length > 0}
        <div class="flex gap-1.5">
          <button
            onclick={() => navigator.clipboard.writeText(displayLog.join("\n"))}
            class="rounded border text-[0.62rem] px-1.5 py-px cursor-pointer transition-colors
                   border-slate-300 text-slate-400 hover:border-slate-400 hover:text-slate-600 bg-transparent
                   dark:border-[#333] dark:text-muted dark:hover:border-[#666] dark:hover:text-pale">
            Copy
          </button>
          <button
            onclick={() => { if (logJobId !== null) jobLogs.set(logJobId, []); }}
            class="rounded border text-[0.62rem] px-1.5 py-px cursor-pointer transition-colors
                   border-slate-300 text-slate-400 hover:border-slate-400 hover:text-slate-600 bg-transparent
                   dark:border-[#333] dark:text-muted dark:hover:border-[#666] dark:hover:text-pale">
            Clear
          </button>
        </div>
      {/if}
    </div>
    {#if logJobId !== null && runningJobIds.has(logJobId)}
      {@const bars = jobProgress.get(logJobId) ?? []}
      {#if bars.length > 0}
        <div class="flex flex-col gap-1.5 pb-2 mb-2 border-b border-slate-200 dark:border-well shrink-0">
          {#each bars as bar}
            {@const pct = bar.total > 0 ? bar.current / bar.total : 0}
            <div class="flex items-center gap-2 text-[0.65rem]">
              <span class="w-16 shrink-0 text-slate-400 dark:text-muted capitalize">{bar.stage}</span>
              <div class="flex-1 h-1.5 rounded-full bg-slate-200 dark:bg-navy overflow-hidden">
                <div class="h-full rounded-full bg-brand transition-all duration-300"
                     style="width:{(pct * 100).toFixed(1)}%">
                </div>
              </div>
              <span class="w-16 shrink-0 text-right tabular-nums text-slate-400 dark:text-muted">
                {bar.current}/{bar.total}
              </span>
            </div>
          {/each}
        </div>
      {/if}
    {/if}
    <pre
      class="flex-1 overflow-y-auto m-0 text-[0.7rem] leading-relaxed whitespace-pre-wrap break-all font-mono
             text-slate-600 dark:text-[#b0c4de]"
      bind:this={logEl}>{displayLog.join("\n")}</pre>
  </div>

</div>

<!-- Context menu -->
{#if ctxMenu}
  <!-- svelte-ignore a11y_click_events_have_key_events a11y_no_static_element_interactions -->
  <div class="fixed inset-0 z-[99]" onclick={closeCtxMenu}></div>
  <div class="fixed z-[100] rounded-xl py-1 min-w-[150px] shadow-2xl border
              bg-white border-slate-200
              dark:bg-[#1a2744] dark:border-[#2a4a8a]"
       style="left:{ctxMenu.x}px;top:{ctxMenu.y}px">
    <button
      onclick={() => startRename(ctxMenu!.job)}
      class="block w-full text-left border-none px-3.5 py-1.5 text-[0.82rem] cursor-pointer transition-colors
             bg-transparent text-slate-700 hover:bg-slate-100
             dark:text-pale dark:hover:bg-[#2a4a8a]">
      Rename
    </button>
    <button
      onclick={() => reloadJob(ctxMenu!.job)}
      class="block w-full text-left border-none px-3.5 py-1.5 text-[0.82rem] cursor-pointer transition-colors
             bg-transparent text-slate-700 hover:bg-slate-100
             dark:text-pale dark:hover:bg-[#2a4a8a]">
      Reload settings
    </button>
    <button
      onclick={() => deleteJob(ctxMenu!.job)}
      class="block w-full text-left border-none px-3.5 py-1.5 text-[0.82rem] cursor-pointer transition-colors
             bg-transparent text-red-500 hover:bg-red-50
             dark:text-brand dark:hover:bg-[#3a1a1a]">
      Delete from history
    </button>
  </div>
{/if}
