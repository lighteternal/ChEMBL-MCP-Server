#!/usr/bin/env node
import { Server } from '@modelcontextprotocol/sdk/server/index.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import {
  CallToolRequestSchema,
  ErrorCode,
  ListResourcesRequestSchema,
  ListResourceTemplatesRequestSchema,
  ListToolsRequestSchema,
  McpError,
  ReadResourceRequestSchema,
} from '@modelcontextprotocol/sdk/types.js';
import axios, { AxiosInstance } from 'axios';

// ChEMBL API interfaces
interface CompoundSearchResult {
  molecule_chembl_id: string;
  pref_name?: string;
  molecule_type: string;
  molecule_structures: {
    canonical_smiles?: string;
    standard_inchi?: string;
    standard_inchi_key?: string;
  };
  molecule_properties: {
    molecular_weight?: number;
    alogp?: number;
    hbd?: number;
    hba?: number;
    psa?: number;
    rtb?: number;
    ro3_pass?: string;
    num_ro5_violations?: number;
  };
}

interface TargetInfo {
  target_chembl_id: string;
  pref_name: string;
  target_type: string;
  organism: string;
  species_group_flag: boolean;
  target_components?: Array<{
    component_id: number;
    component_type: string;
    accession?: string;
    sequence?: string;
  }>;
}

interface ActivityData {
  activity_id: number;
  assay_chembl_id: string;
  molecule_chembl_id: string;
  target_chembl_id: string;
  standard_type?: string;
  standard_value?: number;
  standard_units?: string;
  standard_relation?: string;
  activity_comment?: string;
}

interface AssayInfo {
  assay_chembl_id: string;
  description: string;
  assay_type: string;
  assay_organism?: string;
  assay_strain?: string;
  assay_tissue?: string;
  assay_cell_type?: string;
  assay_subcellular_fraction?: string;
  target_chembl_id?: string;
  confidence_score?: number;
}

// Type guards and validation functions
const isValidCompoundSearchArgs = (
  args: any
): args is { query: string; limit?: number; offset?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.query === 'string' &&
    args.query.length > 0 &&
    (args.limit === undefined || (typeof args.limit === 'number' && args.limit > 0 && args.limit <= 1000)) &&
    (args.offset === undefined || (typeof args.offset === 'number' && args.offset >= 0))
  );
};

const isValidChemblIdArgs = (
  args: any
): args is { chembl_id: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.chembl_id === 'string' &&
    args.chembl_id.length > 0
  );
};

const isValidSimilaritySearchArgs = (
  args: any
): args is { smiles: string; similarity?: number; limit?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.smiles === 'string' &&
    args.smiles.length > 0 &&
    (args.similarity === undefined || (typeof args.similarity === 'number' && args.similarity >= 0 && args.similarity <= 1)) &&
    (args.limit === undefined || (typeof args.limit === 'number' && args.limit > 0 && args.limit <= 1000))
  );
};

const isValidSubstructureSearchArgs = (
  args: any
): args is { smiles: string; limit?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.smiles === 'string' &&
    args.smiles.length > 0 &&
    (args.limit === undefined || (typeof args.limit === 'number' && args.limit > 0 && args.limit <= 1000))
  );
};

const isValidActivitySearchArgs = (
  args: any
): args is { target_chembl_id?: string; assay_chembl_id?: string; molecule_chembl_id?: string; activity_type?: string; limit?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (args.target_chembl_id === undefined || typeof args.target_chembl_id === 'string') &&
    (args.assay_chembl_id === undefined || typeof args.assay_chembl_id === 'string') &&
    (args.molecule_chembl_id === undefined || typeof args.molecule_chembl_id === 'string') &&
    (args.activity_type === undefined || typeof args.activity_type === 'string') &&
    (args.limit === undefined || (typeof args.limit === 'number' && args.limit > 0 && args.limit <= 1000)) &&
    (args.target_chembl_id !== undefined || args.assay_chembl_id !== undefined || args.molecule_chembl_id !== undefined)
  );
};

const isValidPropertyFilterArgs = (
  args: any
): args is {
    min_mw?: number;
    max_mw?: number;
    min_logp?: number;
    max_logp?: number;
    max_hbd?: number;
    max_hba?: number;
    limit?: number
  } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (args.min_mw === undefined || (typeof args.min_mw === 'number' && args.min_mw >= 0)) &&
    (args.max_mw === undefined || (typeof args.max_mw === 'number' && args.max_mw >= 0)) &&
    (args.min_logp === undefined || typeof args.min_logp === 'number') &&
    (args.max_logp === undefined || typeof args.max_logp === 'number') &&
    (args.max_hbd === undefined || (typeof args.max_hbd === 'number' && args.max_hbd >= 0)) &&
    (args.max_hba === undefined || (typeof args.max_hba === 'number' && args.max_hba >= 0)) &&
    (args.limit === undefined || (typeof args.limit === 'number' && args.limit > 0 && args.limit <= 1000))
  );
};

const isValidBatchArgs = (
  args: any
): args is { chembl_ids: string[] } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    Array.isArray(args.chembl_ids) &&
    args.chembl_ids.length > 0 &&
    args.chembl_ids.length <= 50 &&
    args.chembl_ids.every((id: any) => typeof id === 'string' && id.length > 0)
  );
};

const isValidInchiSearchArgs = (
  args: any
): args is { inchi: string; limit?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.inchi === 'string' &&
    args.inchi.length > 0 &&
    (args.limit === undefined || (typeof args.limit === 'number' && args.limit > 0 && args.limit <= 1000))
  );
};

class ChEMBLServer {
  private server: Server;
  private apiClient: AxiosInstance;

  constructor() {
    this.server = new Server(
      {
        name: 'chembl-server',
        version: '1.0.0',
      },
      {
        capabilities: {
          resources: {},
          tools: {},
        },
      }
    );

    // Initialize ChEMBL API client
    this.apiClient = axios.create({
      baseURL: 'https://www.ebi.ac.uk/chembl/api/data',
      timeout: 30000,
      headers: {
        'User-Agent': 'ChEMBL-MCP-Server/1.0.0',
        'Accept': 'application/json',
      },
    });

    this.setupResourceHandlers();
    this.setupToolHandlers();

    // Error handling
    this.server.onerror = (error: any) => console.error('[MCP Error]', error);
    process.on('SIGINT', async () => {
      await this.server.close();
      process.exit(0);
    });
  }

  private setupResourceHandlers() {
    // List available resource templates
    this.server.setRequestHandler(
      ListResourceTemplatesRequestSchema,
      async () => ({
        resourceTemplates: [
          {
            uriTemplate: 'chembl://compound/{chembl_id}',
            name: 'ChEMBL compound entry',
            mimeType: 'application/json',
            description: 'Complete compound information for a ChEMBL ID',
          },
          {
            uriTemplate: 'chembl://target/{chembl_id}',
            name: 'ChEMBL target entry',
            mimeType: 'application/json',
            description: 'Complete target information for a ChEMBL target ID',
          },
          {
            uriTemplate: 'chembl://assay/{chembl_id}',
            name: 'ChEMBL assay entry',
            mimeType: 'application/json',
            description: 'Complete assay information for a ChEMBL assay ID',
          },
          {
            uriTemplate: 'chembl://activity/{activity_id}',
            name: 'ChEMBL activity entry',
            mimeType: 'application/json',
            description: 'Bioactivity measurement data for an activity ID',
          },
          {
            uriTemplate: 'chembl://search/{query}',
            name: 'ChEMBL search results',
            mimeType: 'application/json',
            description: 'Search results for compounds matching the query',
          },
        ],
      })
    );

    // Handle resource requests
    this.server.setRequestHandler(
      ReadResourceRequestSchema,
      async (request: any) => {
        const uri = request.params.uri;

        // Handle compound info requests
        const compoundMatch = uri.match(/^chembl:\/\/compound\/([A-Z0-9]+)$/);
        if (compoundMatch) {
          const chemblId = compoundMatch[1];
          try {
            const response = await this.apiClient.get(`/molecule/${chemblId}.json`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch compound ${chemblId}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle target info requests
        const targetMatch = uri.match(/^chembl:\/\/target\/([A-Z0-9]+)$/);
        if (targetMatch) {
          const chemblId = targetMatch[1];
          try {
            const response = await this.apiClient.get(`/target/${chemblId}.json`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch target ${chemblId}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle assay info requests
        const assayMatch = uri.match(/^chembl:\/\/assay\/([A-Z0-9]+)$/);
        if (assayMatch) {
          const chemblId = assayMatch[1];
          try {
            const response = await this.apiClient.get(`/assay/${chemblId}.json`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch assay ${chemblId}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle activity info requests
        const activityMatch = uri.match(/^chembl:\/\/activity\/([0-9]+)$/);
        if (activityMatch) {
          const activityId = activityMatch[1];
          try {
            const response = await this.apiClient.get(`/activity/${activityId}.json`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch activity ${activityId}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle search requests
        const searchMatch = uri.match(/^chembl:\/\/search\/(.+)$/);
        if (searchMatch) {
          const query = decodeURIComponent(searchMatch[1]);
          try {
            const response = await this.apiClient.get('/molecule/search.json', {
              params: {
                q: query,
                limit: 25,
              },
            });

            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to search compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        throw new McpError(
          ErrorCode.InvalidRequest,
          `Invalid URI format: ${uri}`
        );
      }
    );
  }

  private setupToolHandlers() {
    this.server.setRequestHandler(ListToolsRequestSchema, async () => ({
      tools: [
        // Core Chemical Search & Retrieval (5 tools)
        {
          name: 'search_compounds',
          description: 'Search ChEMBL database for compounds by name, synonym, or identifier',
          inputSchema: {
            type: 'object',
            properties: {
              query: { type: 'string', description: 'Search query (compound name, synonym, or identifier)' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
              offset: { type: 'number', description: 'Number of results to skip (default: 0)', minimum: 0 },
            },
            required: ['query'],
          },
        },
        {
          name: 'get_compound_info',
          description: 'Get detailed information for a specific compound by ChEMBL ID',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID (e.g., CHEMBL59)' },
            },
            required: ['chembl_id'],
          },
        },
        {
          name: 'search_by_inchi',
          description: 'Search for compounds by InChI key or InChI string',
          inputSchema: {
            type: 'object',
            properties: {
              inchi: { type: 'string', description: 'InChI key or InChI string' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['inchi'],
          },
        },
        {
          name: 'get_compound_structure',
          description: 'Retrieve chemical structure information in various formats',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
              format: { type: 'string', enum: ['smiles', 'inchi', 'molfile', 'sdf'], description: 'Structure format (default: smiles)' },
            },
            required: ['chembl_id'],
          },
        },
        {
          name: 'search_similar_compounds',
          description: 'Find chemically similar compounds using Tanimoto similarity',
          inputSchema: {
            type: 'object',
            properties: {
              smiles: { type: 'string', description: 'SMILES string of the query molecule' },
              similarity: { type: 'number', description: 'Similarity threshold (0-1, default: 0.7)', minimum: 0, maximum: 1 },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['smiles'],
          },
        },
        // Target Analysis & Drug Discovery (5 tools)
        {
          name: 'search_targets',
          description: 'Search for biological targets by name or type',
          inputSchema: {
            type: 'object',
            properties: {
              query: { type: 'string', description: 'Target name or search query' },
              target_type: { type: 'string', description: 'Target type filter (e.g., SINGLE PROTEIN, PROTEIN COMPLEX)' },
              organism: { type: 'string', description: 'Organism filter' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['query'],
          },
        },
        {
          name: 'get_target_info',
          description: 'Get detailed information for a specific target by ChEMBL target ID',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL target ID (e.g., CHEMBL2095173)' },
            },
            required: ['chembl_id'],
          },
        },
        {
          name: 'get_target_compounds',
          description: 'Get compounds tested against a specific target',
          inputSchema: {
            type: 'object',
            properties: {
              target_chembl_id: { type: 'string', description: 'ChEMBL target ID' },
              activity_type: { type: 'string', description: 'Activity type filter (e.g., IC50, Ki, Kd)' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['target_chembl_id'],
          },
        },
        {
          name: 'search_by_uniprot',
          description: 'Find ChEMBL targets by UniProt accession',
          inputSchema: {
            type: 'object',
            properties: {
              uniprot_id: { type: 'string', description: 'UniProt accession number' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['uniprot_id'],
          },
        },
        {
          name: 'get_target_pathways',
          description: 'Get biological pathways associated with a target',
          inputSchema: {
            type: 'object',
            properties: {
              target_chembl_id: { type: 'string', description: 'ChEMBL target ID' },
            },
            required: ['target_chembl_id'],
          },
        },
        // Bioactivity & Assay Data (5 tools)
        {
          name: 'search_activities',
          description: 'Search bioactivity measurements and assay results',
          inputSchema: {
            type: 'object',
            properties: {
              target_chembl_id: { type: 'string', description: 'ChEMBL target ID filter' },
              assay_chembl_id: { type: 'string', description: 'ChEMBL assay ID filter' },
              molecule_chembl_id: { type: 'string', description: 'ChEMBL compound ID filter' },
              activity_type: { type: 'string', description: 'Activity type (e.g., IC50, Ki, EC50)' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: [],
          },
        },
        {
          name: 'get_assay_info',
          description: 'Get detailed information for a specific assay by ChEMBL assay ID',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL assay ID (e.g., CHEMBL1217643)' },
            },
            required: ['chembl_id'],
          },
        },
        {
          name: 'search_by_activity_type',
          description: 'Find bioactivity data by specific activity type and value range',
          inputSchema: {
            type: 'object',
            properties: {
              activity_type: { type: 'string', description: 'Activity type (e.g., IC50, Ki, EC50, Kd)' },
              min_value: { type: 'number', description: 'Minimum activity value' },
              max_value: { type: 'number', description: 'Maximum activity value' },
              units: { type: 'string', description: 'Units filter (e.g., nM, uM)' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['activity_type'],
          },
        },
        {
          name: 'get_dose_response',
          description: 'Get dose-response data and activity profiles for compounds',
          inputSchema: {
            type: 'object',
            properties: {
              molecule_chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
              target_chembl_id: { type: 'string', description: 'ChEMBL target ID (optional filter)' },
            },
            required: ['molecule_chembl_id'],
          },
        },
        {
          name: 'compare_activities',
          description: 'Compare bioactivity data across multiple compounds or targets',
          inputSchema: {
            type: 'object',
            properties: {
              molecule_chembl_ids: { type: 'array', items: { type: 'string' }, description: 'Array of ChEMBL compound IDs (2-10)', minItems: 2, maxItems: 10 },
              target_chembl_id: { type: 'string', description: 'ChEMBL target ID for comparison' },
              activity_type: { type: 'string', description: 'Activity type for comparison' },
            },
            required: ['molecule_chembl_ids'],
          },
        },
        // Drug Development & Clinical Data (4 tools)
        {
          name: 'search_drugs',
          description: 'Search for approved drugs and clinical candidates',
          inputSchema: {
            type: 'object',
            properties: {
              query: { type: 'string', description: 'Drug name or search query' },
              development_phase: { type: 'string', description: 'Development phase filter (e.g., Approved, Phase III)' },
              therapeutic_area: { type: 'string', description: 'Therapeutic area filter' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['query'],
          },
        },
        {
          name: 'get_drug_info',
          description: 'Get drug development status and clinical trial information',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
            },
            required: ['chembl_id'],
          },
        },
        {
          name: 'search_drug_indications',
          description: 'Search for therapeutic indications and disease areas',
          inputSchema: {
            type: 'object',
            properties: {
              indication: { type: 'string', description: 'Disease or indication search term' },
              drug_type: { type: 'string', description: 'Drug type filter (e.g., Small molecule, Antibody)' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['indication'],
          },
        },
        {
          name: 'get_mechanism_of_action',
          description: 'Get mechanism of action and target interaction data',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
            },
            required: ['chembl_id'],
          },
        },
        // Chemical Property Analysis (4 tools)
        {
          name: 'analyze_admet_properties',
          description: 'Analyze ADMET properties (Absorption, Distribution, Metabolism, Excretion, Toxicity)',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
            },
            required: ['chembl_id'],
          },
        },
        {
          name: 'calculate_descriptors',
          description: 'Calculate molecular descriptors and physicochemical properties',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
              smiles: { type: 'string', description: 'SMILES string (alternative to ChEMBL ID)' },
            },
            required: [],
          },
        },
        {
          name: 'predict_solubility',
          description: 'Predict aqueous solubility and permeability properties',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
              smiles: { type: 'string', description: 'SMILES string (alternative to ChEMBL ID)' },
            },
            required: [],
          },
        },
        {
          name: 'assess_drug_likeness',
          description: 'Assess drug-likeness using Lipinski Rule of Five and other metrics',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound ID' },
              smiles: { type: 'string', description: 'SMILES string (alternative to ChEMBL ID)' },
            },
            required: [],
          },
        },
        // Advanced Search & Cross-References (4 tools)
        {
          name: 'substructure_search',
          description: 'Find compounds containing specific substructures',
          inputSchema: {
            type: 'object',
            properties: {
              smiles: { type: 'string', description: 'SMILES string of the substructure query' },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: ['smiles'],
          },
        },
        {
          name: 'batch_compound_lookup',
          description: 'Process multiple ChEMBL IDs efficiently',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_ids: { type: 'array', items: { type: 'string' }, description: 'Array of ChEMBL compound IDs (1-50)', minItems: 1, maxItems: 50 },
            },
            required: ['chembl_ids'],
          },
        },
        {
          name: 'get_external_references',
          description: 'Get links to external databases (PubChem, DrugBank, PDB, etc.)',
          inputSchema: {
            type: 'object',
            properties: {
              chembl_id: { type: 'string', description: 'ChEMBL compound or target ID' },
            },
            required: ['chembl_id'],
          },
        },
        {
          name: 'advanced_search',
          description: 'Complex queries with multiple chemical and biological filters',
          inputSchema: {
            type: 'object',
            properties: {
              min_mw: { type: 'number', description: 'Minimum molecular weight (Da)', minimum: 0 },
              max_mw: { type: 'number', description: 'Maximum molecular weight (Da)', minimum: 0 },
              min_logp: { type: 'number', description: 'Minimum LogP value' },
              max_logp: { type: 'number', description: 'Maximum LogP value' },
              max_hbd: { type: 'number', description: 'Maximum hydrogen bond donors', minimum: 0 },
              max_hba: { type: 'number', description: 'Maximum hydrogen bond acceptors', minimum: 0 },
              limit: { type: 'number', description: 'Number of results to return (1-1000, default: 25)', minimum: 1, maximum: 1000 },
            },
            required: [],
          },
        },
      ],
    }));

    this.server.setRequestHandler(CallToolRequestSchema, async (request: any) => {
      const { name, arguments: args } = request.params;

      try {
        switch (name) {
          // Core Chemical Search & Retrieval
          case 'search_compounds':
            return await this.handleSearchCompounds(args);
          case 'get_compound_info':
            return await this.handleGetCompoundInfo(args);
          case 'search_by_inchi':
            return await this.handleSearchByInchi(args);
          case 'get_compound_structure':
            return await this.handleGetCompoundStructure(args);
          case 'search_similar_compounds':
            return await this.handleSearchSimilarCompounds(args);
          // Target Analysis & Drug Discovery
          case 'search_targets':
            return await this.handleSearchTargets(args);
          case 'get_target_info':
            return await this.handleGetTargetInfo(args);
          case 'get_target_compounds':
            return await this.handleGetTargetCompounds(args);
          case 'search_by_uniprot':
            return await this.handleSearchByUniprot(args);
          case 'get_target_pathways':
            return await this.handleGetTargetPathways(args);
          // Bioactivity & Assay Data
          case 'search_activities':
            return await this.handleSearchActivities(args);
          case 'get_assay_info':
            return await this.handleGetAssayInfo(args);
          case 'search_by_activity_type':
            return await this.handleSearchByActivityType(args);
          case 'get_dose_response':
            return await this.handleGetDoseResponse(args);
          case 'compare_activities':
            return await this.handleCompareActivities(args);
          // Drug Development & Clinical Data
          case 'search_drugs':
            return await this.handleSearchDrugs(args);
          case 'get_drug_info':
            return await this.handleGetDrugInfo(args);
          case 'search_drug_indications':
            return await this.handleSearchDrugIndications(args);
          case 'get_mechanism_of_action':
            return await this.handleGetMechanismOfAction(args);
          // Chemical Property Analysis
          case 'analyze_admet_properties':
            return await this.handleAnalyzeAdmetProperties(args);
          case 'calculate_descriptors':
            return await this.handleCalculateDescriptors(args);
          case 'predict_solubility':
            return await this.handlePredictSolubility(args);
          case 'assess_drug_likeness':
            return await this.handleAssessDrugLikeness(args);
          // Advanced Search & Cross-References
          case 'substructure_search':
            return await this.handleSubstructureSearch(args);
          case 'batch_compound_lookup':
            return await this.handleBatchCompoundLookup(args);
          case 'get_external_references':
            return await this.handleGetExternalReferences(args);
          case 'advanced_search':
            return await this.handleAdvancedSearch(args);
          default:
            throw new McpError(
              ErrorCode.MethodNotFound,
              `Unknown tool: ${name}`
            );
        }
      } catch (error) {
        return {
          content: [
            {
              type: 'text',
              text: `Error executing tool ${name}: ${error instanceof Error ? error.message : 'Unknown error'}`,
            },
          ],
          isError: true,
        };
      }
    });
  }

  // Core Chemical Search & Retrieval handlers
  private async handleSearchCompounds(args: any) {
    if (!isValidCompoundSearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid compound search arguments');
    }

    try {
      const response = await this.apiClient.get('/molecule/search.json', {
        params: {
          q: args.query,
          limit: args.limit || 25,
          offset: args.offset || 0,
        },
      });

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetCompoundInfo(args: any) {
    if (!isValidChemblIdArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid ChEMBL ID arguments');
    }

    try {
      const response = await this.apiClient.get(`/molecule/${args.chembl_id}.json`);
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get compound info: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // Simplified placeholder implementations for the remaining tools
  private async handleSearchByInchi(args: any) {
    if (!isValidInchiSearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments for InChI search.');
    }

    const params: any = {
        limit: args.limit || 25,
    };

    if (args.inchi.startsWith('InChI=')) {
        params['molecule_structures__standard_inchi'] = args.inchi;
    } else {
        params['molecule_structures__standard_inchi_key'] = args.inchi;
    }

    try {
      const response = await this.apiClient.get('/molecule.json', { params });
      return {
        content: [{ type: 'text', text: JSON.stringify(response.data, null, 2) }],
      };
    } catch (error: any) {
      throw new McpError(ErrorCode.InternalError, `ChEMBL API error: ${error.message}`);
    }
  }

  private async handleGetCompoundStructure(args: any) {
    if (!args || typeof args.chembl_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments');
    }

    try {
      const response = await this.apiClient.get(`/molecule/${args.chembl_id}.json`);
      const compound = response.data;

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({
              chembl_id: compound.molecule_chembl_id,
              structures: compound.molecule_structures || {},
              requested_format: args.format || 'smiles'
            }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Failed to get structure: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  private async handleSearchSimilarCompounds(args: any) {
    if (!isValidSimilaritySearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid similarity search arguments');
    }

    const similarity = Math.round((args.similarity || 0.7) * 100);
    const smiles = encodeURIComponent(args.smiles);

    try {
      const response = await this.apiClient.get(`/similarity/${smiles}/${similarity}.json`, {
        params: {
          limit: args.limit || 25,
        },
      });

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search similar compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSearchTargets(args: any) {
    try {
      const response = await this.apiClient.get('/target/search.json', {
        params: { q: args.query, limit: args.limit || 25 },
      });
      return { content: [{ type: 'text', text: JSON.stringify(response.data, null, 2) }] };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Failed to search targets: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  private async handleGetTargetInfo(args: any) {
    if (!isValidChemblIdArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments');
    }

    try {
      const response = await this.apiClient.get(`/target/${args.chembl_id}.json`);
      return { content: [{ type: 'text', text: JSON.stringify(response.data, null, 2) }] };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Failed to get target info: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  // Placeholder implementations for remaining tools
  private async handleGetTargetCompounds(args: any) {
    if (!args || typeof args.target_chembl_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: target_chembl_id is required');
    }

    try {
      const params: any = {
        target_chembl_id: args.target_chembl_id,
        limit: args.limit || 25,
      };

      if (args.activity_type) {
        params.standard_type = args.activity_type;
      }

      const response = await this.apiClient.get('/activity.json', { params });
      
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get target compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSearchByUniprot(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'UniProt search not yet implemented', args }, null, 2) }] };
  }

  private async handleGetTargetPathways(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Target pathways not yet implemented', args }, null, 2) }] };
  }

  private async handleSearchActivities(args: any) {
    try {
      const params: any = { limit: args.limit || 25 };
      if (args.target_chembl_id) params.target_chembl_id = args.target_chembl_id;
      if (args.molecule_chembl_id) params.molecule_chembl_id = args.molecule_chembl_id;
      if (args.activity_type) params.standard_type = args.activity_type;

      const response = await this.apiClient.get('/activity.json', { params });
      return { content: [{ type: 'text', text: JSON.stringify(response.data, null, 2) }] };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Failed to search activities: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  private async handleGetAssayInfo(args: any) {
    if (!isValidChemblIdArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments');
    }

    try {
      const response = await this.apiClient.get(`/assay/${args.chembl_id}.json`);
      return { content: [{ type: 'text', text: JSON.stringify(response.data, null, 2) }] };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Failed to get assay info: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  // Remaining placeholder implementations
  private async handleSearchByActivityType(args: any) {
    if (!args || typeof args.activity_type !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: activity_type is required');
    }

    try {
      const params: any = {
        standard_type: args.activity_type,
        limit: Math.min(args.limit || 10, 20), // Limit to max 20 results
      };

      if (args.min_value !== undefined) {
        params.standard_value__gte = args.min_value;
      }
      if (args.max_value !== undefined) {
        params.standard_value__lte = args.max_value;
      }
      if (args.units) {
        params.standard_units = args.units;
      }

      const response = await this.apiClient.get('/activity.json', { params });
      
      // Extract and clean essential data only
      const activities = response.data.activities || [];
      const cleanedActivities = activities.slice(0, 10).map((activity: any) => ({
        activity_id: activity.activity_id,
        molecule_chembl_id: activity.molecule_chembl_id,
        target_chembl_id: activity.target_chembl_id,
        assay_chembl_id: activity.assay_chembl_id,
        standard_type: activity.standard_type,
        standard_value: activity.standard_value,
        standard_units: activity.standard_units,
        pchembl_value: activity.pchembl_value,
        assay_description: activity.assay_description?.substring(0, 200) || '', // Truncate long descriptions
        document_journal: activity.document_journal,
        document_year: activity.document_year,
      }));

      const summary = {
        total_count: response.data.page_meta?.total_count || activities.length,
        returned_count: cleanedActivities.length,
        activity_type: args.activity_type,
        activities: cleanedActivities,
      };
      
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(summary, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search by activity type: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetDoseResponse(args: any) {
    if (!args || typeof args.molecule_chembl_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: molecule_chembl_id is required');
    }

    try {
      const params: any = {
        molecule_chembl_id: args.molecule_chembl_id,
        limit: args.limit || 25,
      };

      if (args.target_chembl_id) {
        params.target_chembl_id = args.target_chembl_id;
      }

      // Focus on dose-response related activity types
      params.standard_type__in = 'IC50,EC50,Ki,Kd,GI50,LC50,LD50';

      const response = await this.apiClient.get('/activity.json', { params });
      
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get dose response data: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleCompareActivities(args: any) {
    if (!args || !Array.isArray(args.molecule_chembl_ids) || args.molecule_chembl_ids.length < 2) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: molecule_chembl_ids array with at least 2 compounds is required');
    }

    try {
      const results = [];
      
      for (const chemblId of args.molecule_chembl_ids.slice(0, 5)) { // Limit to 5 compounds
        const params: any = {
          molecule_chembl_id: chemblId,
          limit: 10, // Reduced limit per compound
        };

        if (args.target_chembl_id) {
          params.target_chembl_id = args.target_chembl_id;
        }
        if (args.activity_type) {
          params.standard_type = args.activity_type;
        }

        try {
          const response = await this.apiClient.get('/activity.json', { params });
          const activities = response.data.activities || [];
          
          // Extract only essential activity data
          const cleanedActivities = activities.slice(0, 5).map((activity: any) => ({
            activity_id: activity.activity_id,
            target_chembl_id: activity.target_chembl_id,
            assay_chembl_id: activity.assay_chembl_id,
            standard_type: activity.standard_type,
            standard_value: activity.standard_value,
            standard_units: activity.standard_units,
            pchembl_value: activity.pchembl_value,
            assay_description: activity.assay_description?.substring(0, 100) || '',
          }));

          results.push({
            molecule_chembl_id: chemblId,
            activity_count: activities.length,
            top_activities: cleanedActivities,
            success: true
          });
        } catch (error) {
          results.push({
            molecule_chembl_id: chemblId,
            error: error instanceof Error ? error.message : 'Unknown error',
            success: false
          });
        }
      }

      const summary = {
        comparison_type: args.activity_type || 'all_activities',
        target_filter: args.target_chembl_id || 'all_targets',
        compounds_compared: results.length,
        comparison_results: results,
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(summary, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to compare activities: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSearchDrugs(args: any) {
    if (!isValidCompoundSearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments for drug search.');
    }
    try {
      const response = await this.apiClient.get('/drug.json', {
        params: {
          pref_name__icontains: args.query,
          limit: args.limit || 10,
          offset: args.offset || 0,
        },
      });
      return {
        content: [{ type: 'text', text: JSON.stringify(response.data, null, 2) }],
      };
    } catch (error: any) {
      throw new McpError(ErrorCode.InternalError, `ChEMBL API error: ${error.message}`);
    }
  }

  private async handleGetDrugInfo(args: any) {
    if (!isValidChemblIdArgs(args)) {
        throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments for getting drug info.');
    }
    return this.handleGetCompoundInfo(args);
  }

  private async handleSearchDrugIndications(args: any) {
    if (typeof args.indication !== 'string' || args.indication.length === 0) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments for drug indication search.');
    }
    try {
      const response = await this.apiClient.get('/drug_indication.json', {
        params: {
          efo_term__icontains: args.indication,
          limit: args.limit || 10,
          offset: args.offset || 0,
        },
      });
      return {
        content: [{ type: 'text', text: JSON.stringify(response.data, null, 2) }],
      };
    } catch (error: any) {
      throw new McpError(ErrorCode.InternalError, `ChEMBL API error: ${error.message}`);
    }
  }

  private async handleGetMechanismOfAction(args: any) {
    if (!args || typeof args.chembl_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: chembl_id is required');
    }

    try {
      const response = await this.apiClient.get('/mechanism.json', {
        params: {
          molecule_chembl_id: args.chembl_id,
          limit: 25,
        },
      });
      
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get mechanism of action: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleAnalyzeAdmetProperties(args: any) {
    if (!args || typeof args.chembl_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: chembl_id is required');
    }

    try {
      const response = await this.apiClient.get(`/molecule/${args.chembl_id}.json`);
      const compound = response.data;

      // Extract ADMET-relevant properties
      const admetAnalysis = {
        chembl_id: compound.molecule_chembl_id,
        pref_name: compound.pref_name,
        molecular_properties: compound.molecule_properties,
        admet_assessment: {
          absorption: {
            molecular_weight: compound.molecule_properties?.molecular_weight,
            lipophilicity_alogp: compound.molecule_properties?.alogp,
            polar_surface_area: compound.molecule_properties?.psa,
            hbd_count: compound.molecule_properties?.hbd,
            hba_count: compound.molecule_properties?.hba,
            rotatable_bonds: compound.molecule_properties?.rtb,
          },
          drug_likeness: {
            lipinski_violations: compound.molecule_properties?.num_ro5_violations,
            ro3_pass: compound.molecule_properties?.ro3_pass,
          },
          structural_alerts: compound.compound_structural_alerts || [],
        },
        atc_classifications: compound.atc_classifications || [],
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(admetAnalysis, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to analyze ADMET properties: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleCalculateDescriptors(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Descriptor calculation not yet implemented', args }, null, 2) }] };
  }

  private async handlePredictSolubility(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Solubility prediction not yet implemented', args }, null, 2) }] };
  }

  private async handleAssessDrugLikeness(args: any) {
    if (!args || (typeof args.chembl_id !== 'string' && typeof args.smiles !== 'string')) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: chembl_id or smiles is required');
    }

    try {
      let compound;
      
      if (args.chembl_id) {
        const response = await this.apiClient.get(`/molecule/${args.chembl_id}.json`);
        compound = response.data;
      } else {
        // For SMILES input, search for the compound first
        const searchResponse = await this.apiClient.get('/molecule.json', {
          params: {
            molecule_structures__canonical_smiles: args.smiles,
            limit: 1,
          },
        });
        
        if (searchResponse.data.molecules && searchResponse.data.molecules.length > 0) {
          compound = searchResponse.data.molecules[0];
        } else {
          throw new Error('Compound not found for the provided SMILES');
        }
      }

      const props = compound.molecule_properties || {};
      
      // Lipinski Rule of Five assessment
      const lipinski = {
        molecular_weight: props.molecular_weight,
        mw_pass: props.molecular_weight <= 500,
        logp: props.alogp,
        logp_pass: props.alogp <= 5,
        hbd: props.hbd,
        hbd_pass: props.hbd <= 5,
        hba: props.hba,
        hba_pass: props.hba <= 10,
        violations: props.num_ro5_violations || 0,
        overall_pass: (props.num_ro5_violations || 0) <= 1,
      };

      // Additional drug-likeness metrics
      const drugLikenessAssessment = {
        chembl_id: compound.molecule_chembl_id,
        pref_name: compound.pref_name,
        smiles: compound.molecule_structures?.canonical_smiles,
        lipinski_rule_of_five: lipinski,
        additional_metrics: {
          polar_surface_area: props.psa,
          psa_pass: props.psa <= 140,
          rotatable_bonds: props.rtb,
          rtb_pass: props.rtb <= 10,
          ro3_pass: props.ro3_pass,
          heavy_atoms: props.heavy_atoms,
          aromatic_rings: props.aromatic_rings,
        },
        overall_drug_likeness: {
          lipinski_compliant: lipinski.overall_pass,
          lead_like: props.ro3_pass === 'Y',
          oral_bioavailability_score: this.calculateOralBioavailabilityScore(props),
        },
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(drugLikenessAssessment, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to assess drug-likeness: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private calculateOralBioavailabilityScore(props: any): string {
    let score = 0;
    if (props.molecular_weight && props.molecular_weight <= 500) score++;
    if (props.alogp && props.alogp <= 5) score++;
    if (props.hbd && props.hbd <= 5) score++;
    if (props.hba && props.hba <= 10) score++;
    if (props.psa && props.psa <= 140) score++;
    if (props.rtb && props.rtb <= 10) score++;
    
    if (score >= 5) return 'High';
    if (score >= 3) return 'Medium';
    return 'Low';
  }

  private async handleSubstructureSearch(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Substructure search not yet implemented', args }, null, 2) }] };
  }

  private async handleBatchCompoundLookup(args: any) {
    if (!isValidBatchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid batch arguments');
    }

    try {
      const results = [];
      for (const chemblId of args.chembl_ids.slice(0, 10)) { // Limit to 10 for demo
        try {
          const response = await this.apiClient.get(`/molecule/${chemblId}.json`);
          results.push({ chembl_id: chemblId, data: response.data, success: true });
        } catch (error) {
          results.push({ chembl_id: chemblId, error: error instanceof Error ? error.message : 'Unknown error', success: false });
        }
      }

      return { content: [{ type: 'text', text: JSON.stringify({ batch_results: results }, null, 2) }] };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Batch lookup failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  private async handleGetExternalReferences(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'External references not yet implemented', args }, null, 2) }] };
  }

  private async handleAdvancedSearch(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Advanced search not yet implemented', args }, null, 2) }] };
  }

  async run() {
    const transport = new StdioServerTransport();
    await this.server.connect(transport);
    console.error('ChEMBL MCP server running on stdio');
  }
}

const server = new ChEMBLServer();
server.run().catch(console.error);
