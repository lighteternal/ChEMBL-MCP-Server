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
    if (!args || typeof args.uniprot_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: uniprot_id is required');
    }

    try {
      // Search for targets that have the specified UniProt accession in their components
      const response = await this.apiClient.get('/target.json', {
        params: {
          target_components__accession: args.uniprot_id,
          limit: args.limit || 25,
        },
      });

      // Extract and format the target data
      const targets = response.data.targets || [];
      const formattedTargets = targets.map((target: any) => ({
        target_chembl_id: target.target_chembl_id,
        pref_name: target.pref_name,
        target_type: target.target_type,
        organism: target.organism,
        species_group_flag: target.species_group_flag,
        target_components: target.target_components?.filter((comp: any) => 
          comp.accession === args.uniprot_id || comp.accession?.includes(args.uniprot_id)
        ) || [],
        cross_references: target.cross_references || [],
      }));

      const result = {
        query: {
          uniprot_id: args.uniprot_id,
          search_type: 'uniprot_accession'
        },
        total_results: response.data.page_meta?.total_count || targets.length,
        results_returned: formattedTargets.length,
        targets: formattedTargets,
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search by UniProt ID: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetTargetPathways(args: any) {
    if (!args || typeof args.target_chembl_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: target_chembl_id is required');
    }

    try {
      // First, get the target information to extract UniProt accessions
      const targetResponse = await this.apiClient.get(`/target/${args.target_chembl_id}.json`);
      const target = targetResponse.data;

      // Extract UniProt accessions from target components
      const uniprotAccessions = target.target_components
        ?.filter((comp: any) => comp.accession && comp.accession.match(/^[A-Z0-9]+$/))
        .map((comp: any) => comp.accession) || [];

             // Get pathway information - ChEMBL doesn't have direct pathway endpoints,
       // but we can get pathway-related cross-references and mechanism data
       const pathwayInfo: any = {
         cross_references: target.cross_references?.filter((ref: any) => 
           ref.xref_src_db && (
             ref.xref_src_db.toLowerCase().includes('reactome') ||
             ref.xref_src_db.toLowerCase().includes('kegg') ||
             ref.xref_src_db.toLowerCase().includes('pathway') ||
             ref.xref_src_db.toLowerCase().includes('biocyc')
           )
         ) || [],
         go_classifications: target.go_classifications || [],
       };

       // Try to get mechanism of action data which may contain pathway information
       try {
         const mechanismResponse = await this.apiClient.get('/mechanism.json', {
           params: {
             target_chembl_id: args.target_chembl_id,
             limit: 50,
           },
         });

         pathwayInfo.mechanisms = mechanismResponse.data.mechanisms?.map((mech: any) => ({
           mechanism_id: mech.mec_id,
           molecule_chembl_id: mech.molecule_chembl_id,
           mechanism_of_action: mech.mechanism_of_action,
           action_type: mech.action_type,
           mechanism_comment: mech.mechanism_comment,
           selectivity_comment: mech.selectivity_comment,
         })) || [];
       } catch (mechanismError) {
         // Mechanism data is optional, continue without it
         pathwayInfo.mechanisms = [];
         pathwayInfo.mechanism_note = 'No mechanism data available';
       }

       const pathwayData = {
         target_chembl_id: target.target_chembl_id,
         pref_name: target.pref_name,
         target_type: target.target_type,
         organism: target.organism,
         uniprot_accessions: uniprotAccessions,
         pathway_information: pathwayInfo,
       };

      // Add a note about pathway data limitations
      const result = {
        ...pathwayData,
        data_note: 'ChEMBL pathway data is limited. For comprehensive pathway information, cross-reference with Reactome, KEGG, or other pathway databases using the provided UniProt accessions.',
        suggested_external_queries: uniprotAccessions.length > 0 ? {
          reactome: `https://reactome.org/content/query?q=${uniprotAccessions[0]}`,
          kegg: `https://www.genome.jp/kegg-bin/search_pathway_text?map=map01100&keyword=${uniprotAccessions[0]}`,
          wikipathways: `https://www.wikipathways.org/index.php/Special:SearchPathways?query=${uniprotAccessions[0]}`,
        } : null,
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get target pathways: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
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
      const structures = compound.molecule_structures || {};

      // Calculate additional descriptors from available data
      const descriptors = {
        chembl_id: compound.molecule_chembl_id,
        pref_name: compound.pref_name,
        smiles: structures.canonical_smiles,
        inchi: structures.standard_inchi,
        inchi_key: structures.standard_inchi_key,
        
        // Basic molecular descriptors
        molecular_properties: {
          molecular_weight: props.molecular_weight,
          exact_mass: props.full_mwt,
          heavy_atoms: props.heavy_atoms,
          num_atoms: props.heavy_atoms ? props.heavy_atoms + (props.hbd || 0) : null,
        },
        
        // Lipophilicity descriptors
        lipophilicity: {
          alogp: props.alogp,
          clogp: props.cx_logp,
          logd: props.cx_logd,
        },
        
        // Hydrogen bonding descriptors
        hydrogen_bonding: {
          hbd: props.hbd, // Hydrogen bond donors
          hba: props.hba, // Hydrogen bond acceptors
          total_hb_sites: (props.hbd || 0) + (props.hba || 0),
        },
        
        // Topological descriptors
        topological: {
          rotatable_bonds: props.rtb,
          aromatic_rings: props.aromatic_rings,
          rings: props.rings,
          aliphatic_rings: props.rings ? props.rings - (props.aromatic_rings || 0) : null,
        },
        
        // Surface area and volume descriptors
        surface_properties: {
          polar_surface_area: props.psa,
          molecular_surface_area: props.molecular_surface_area,
        },
        
        // Drug-likeness related descriptors
        drug_likeness_metrics: {
          lipinski_violations: props.num_ro5_violations,
          ro3_pass: props.ro3_pass,
          bioavailability_score: this.calculateBioavailabilityScore(props),
        },
        
        // Additional calculated properties
        calculated_properties: {
          flexibility_index: props.rtb ? props.rtb / (props.heavy_atoms || 1) : null,
          complexity_index: props.rings && props.aromatic_rings ? 
            (props.rings * 2 + props.aromatic_rings) / (props.heavy_atoms || 1) : null,
          hydrogen_bond_ratio: props.hbd && props.hba ? 
            props.hbd / (props.hbd + props.hba) : null,
        },
        
        // Molecular framework analysis
        framework_analysis: {
          saturation_ratio: props.aromatic_rings && props.rings ? 
            (props.rings - props.aromatic_rings) / props.rings : null,
          ring_density: props.rings && props.heavy_atoms ? 
            props.rings / props.heavy_atoms : null,
        }
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(descriptors, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to calculate descriptors: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private calculateBioavailabilityScore(props: any): number {
    let score = 0;
    let maxScore = 6;
    
    if (props.molecular_weight && props.molecular_weight <= 500) score++;
    if (props.alogp && props.alogp <= 5) score++;
    if (props.hbd && props.hbd <= 5) score++;
    if (props.hba && props.hba <= 10) score++;
    if (props.psa && props.psa <= 140) score++;
    if (props.rtb && props.rtb <= 10) score++;
    
    return score / maxScore;
  }

  private async handlePredictSolubility(args: any) {
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
      const structures = compound.molecule_structures || {};

      // Calculate solubility predictions using established models
      const solubilityPrediction = this.calculateSolubilityPredictions(props);

      const result = {
        chembl_id: compound.molecule_chembl_id,
        pref_name: compound.pref_name,
        smiles: structures.canonical_smiles,
        
        // Input molecular properties
        molecular_properties: {
          molecular_weight: props.molecular_weight,
          alogp: props.alogp,
          hbd: props.hbd,
          hba: props.hba,
          psa: props.psa,
          rotatable_bonds: props.rtb,
        },
        
        // Solubility predictions
        solubility_predictions: {
          ...solubilityPrediction,
          
          // Additional solubility factors
          solubility_factors: {
            lipophilicity_effect: this.assessLipophilicityEffect(props.alogp),
            hydrogen_bonding_effect: this.assessHydrogenBondingEffect(props.hbd, props.hba),
            molecular_size_effect: this.assessMolecularSizeEffect(props.molecular_weight),
            flexibility_effect: this.assessFlexibilityEffect(props.rtb, props.heavy_atoms),
          },
          
          // Overall assessment
          overall_assessment: this.getOverallSolubilityAssessment(solubilityPrediction),
        },
        
        // Recommendations for improving solubility
        recommendations: this.getSolubilityRecommendations(props, solubilityPrediction),
        
        // Methodology note
        methodology: {
          note: 'Predictions based on established QSAR models and molecular descriptors',
          methods_used: [
            'Yalkowsky General Solubility Equation',
            'Lipophilicity-based estimation',
            'Hydrogen bonding contribution',
            'Molecular size correlation'
          ],
          limitations: 'Predictions are estimates and may not account for specific molecular interactions, crystalline forms, or experimental conditions'
        }
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to predict solubility: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private calculateSolubilityPredictions(props: any): any {
    const logP = props.alogp || 0;
    const mw = props.molecular_weight || 0;
    const hbd = props.hbd || 0;
    const hba = props.hba || 0;
    const psa = props.psa || 0;

    // Yalkowsky General Solubility Equation: logS = 0.5 - 0.01(MP - 25) - logP
    // Simplified version without melting point (using average correction)
    const logS_yalkowsky = 0.5 - 0.8 - logP; // -0.8 is average MP correction

    // Ali et al. model: logS = 0.16 - 0.63 * logP + 0.009 * HBD
    const logS_ali = 0.16 - 0.63 * logP + 0.009 * hbd;

    // Delaney (ESOL) simplified: logS = 0.16 - 0.63 * logP - 0.0062 * MW + 0.066 * RB - 0.74
    const rb = props.rtb || 0;
    const logS_delaney = 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74;

    // Average prediction
    const logS_average = (logS_yalkowsky + logS_ali + logS_delaney) / 3;

    // Convert to different units
    const solubility_mgml = Math.pow(10, logS_average) * mw;
    const solubility_molar = Math.pow(10, logS_average);

    return {
      log_solubility: {
        yalkowsky_model: logS_yalkowsky,
        ali_model: logS_ali,
        delaney_model: logS_delaney,
        average: logS_average,
      },
      predicted_solubility: {
        log_s: logS_average,
        molar_solubility: solubility_molar,
        solubility_mg_ml: solubility_mgml,
        solubility_class: this.classifySolubility(logS_average),
      }
    };
  }

  private classifySolubility(logS: number): string {
    if (logS > -1) return 'Very Soluble (>10 mg/mL)';
    if (logS > -2) return 'Soluble (1-10 mg/mL)';
    if (logS > -3) return 'Moderately Soluble (0.1-1 mg/mL)';
    if (logS > -4) return 'Slightly Soluble (0.01-0.1 mg/mL)';
    return 'Poorly Soluble (<0.01 mg/mL)';
  }

  private assessLipophilicityEffect(alogp: number): string {
    if (!alogp) return 'Unknown';
    if (alogp < 0) return 'Hydrophilic - favors solubility';
    if (alogp < 2) return 'Moderate lipophilicity - balanced solubility';
    if (alogp < 4) return 'Lipophilic - reduces aqueous solubility';
    return 'Highly lipophilic - poor aqueous solubility';
  }

  private assessHydrogenBondingEffect(hbd: number, hba: number): string {
    const total = (hbd || 0) + (hba || 0);
    if (total === 0) return 'No hydrogen bonding - poor solubility';
    if (total < 3) return 'Limited hydrogen bonding - moderate effect';
    if (total < 6) return 'Good hydrogen bonding - favors solubility';
    return 'Extensive hydrogen bonding - high solubility potential';
  }

  private assessMolecularSizeEffect(mw: number): string {
    if (!mw) return 'Unknown';
    if (mw < 200) return 'Small molecule - generally more soluble';
    if (mw < 400) return 'Medium size - moderate size effect';
    if (mw < 600) return 'Large molecule - size may limit solubility';
    return 'Very large molecule - significant size hindrance';
  }

  private assessFlexibilityEffect(rtb: number, heavy_atoms: number): string {
    if (!rtb || !heavy_atoms) return 'Unknown';
    const flexibility = rtb / heavy_atoms;
    if (flexibility < 0.1) return 'Rigid structure - may reduce solubility';
    if (flexibility < 0.2) return 'Moderate flexibility - balanced effect';
    return 'Flexible structure - may enhance solubility';
  }

  private getOverallSolubilityAssessment(prediction: any): string {
    const logS = prediction.predicted_solubility.log_s;
    if (logS > -2) return 'GOOD - Expected to have adequate aqueous solubility';
    if (logS > -3) return 'MODERATE - May require formulation optimization';
    if (logS > -4) return 'POOR - Likely to need solubility enhancement strategies';
    return 'VERY POOR - Significant solubility challenges expected';
  }

  private getSolubilityRecommendations(props: any, prediction: any): string[] {
    const recommendations = [];
    const logP = props.alogp || 0;
    const logS = prediction.predicted_solubility.log_s;

    if (logS < -3) {
      recommendations.push('Consider formulation approaches: cyclodextrins, surfactants, or cocrystals');
      recommendations.push('Evaluate salt forms if ionizable groups are present');
      recommendations.push('Consider prodrug strategies to improve solubility');
    }

    if (logP > 3) {
      recommendations.push('High lipophilicity detected - consider adding polar groups');
      recommendations.push('Evaluate hydroxyl, amine, or carboxyl substitutions');
    }

    if ((props.hbd || 0) + (props.hba || 0) < 3) {
      recommendations.push('Limited hydrogen bonding - consider adding H-bond donors/acceptors');
    }

    if (props.molecular_weight > 500) {
      recommendations.push('High molecular weight - consider molecular size reduction');
    }

    if (recommendations.length === 0) {
      recommendations.push('Solubility appears adequate for most applications');
    }

    return recommendations;
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
    if (!isValidSubstructureSearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid substructure search arguments');
    }

    try {
      // ChEMBL supports substructure search via the substructure endpoint
      const smiles = encodeURIComponent(args.smiles);
      
      const response = await this.apiClient.get(`/substructure/${smiles}.json`, {
        params: {
          limit: args.limit || 25,
        },
      });

      // Extract and format the results
      const molecules = response.data.molecules || [];
      const formattedResults = molecules.map((molecule: any) => ({
        molecule_chembl_id: molecule.molecule_chembl_id,
        pref_name: molecule.pref_name,
        molecule_type: molecule.molecule_type,
        max_phase: molecule.max_phase,
        structures: {
          canonical_smiles: molecule.molecule_structures?.canonical_smiles,
          standard_inchi_key: molecule.molecule_structures?.standard_inchi_key,
        },
        properties: {
          molecular_weight: molecule.molecule_properties?.molecular_weight,
          alogp: molecule.molecule_properties?.alogp,
          hbd: molecule.molecule_properties?.hbd,
          hba: molecule.molecule_properties?.hba,
          num_ro5_violations: molecule.molecule_properties?.num_ro5_violations,
        },
        similarity_metrics: {
          contains_query_substructure: true,
          tanimoto_similarity: molecule.similarity || 'N/A',
        }
      }));

      const result = {
        query: {
          query_smiles: args.smiles,
          search_type: 'substructure',
          description: 'Find compounds containing the specified substructure'
        },
        search_results: {
          total_results: response.data.page_meta?.total_count || molecules.length,
          results_returned: formattedResults.length,
          compounds: formattedResults,
        },
        substructure_analysis: {
          query_structure_info: await this.analyzeQueryStructure(args.smiles),
          result_diversity: this.analyzeResultDiversity(formattedResults),
        },
        usage_notes: {
          note: 'Substructure search finds compounds that contain the query structure as a subunit',
          applications: [
            'Scaffold hopping',
            'Lead optimization',
            'Bioisostere identification',
            'Structure-activity relationship analysis'
          ]
        }
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to perform substructure search: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async analyzeQueryStructure(smiles: string): Promise<any> {
    try {
      // Try to get basic information about the query structure
      const searchResponse = await this.apiClient.get('/molecule.json', {
        params: {
          molecule_structures__canonical_smiles: smiles,
          limit: 1,
        },
      });

      if (searchResponse.data.molecules && searchResponse.data.molecules.length > 0) {
        const molecule = searchResponse.data.molecules[0];
        return {
          found_in_chembl: true,
          chembl_id: molecule.molecule_chembl_id,
          pref_name: molecule.pref_name,
          properties: molecule.molecule_properties,
        };
      }
    } catch (error) {
      // Query structure not found in ChEMBL, which is fine
    }

    return {
      found_in_chembl: false,
      smiles: smiles,
      note: 'Query structure not found as exact match in ChEMBL, but used for substructure search',
    };
  }

  private analyzeResultDiversity(results: any[]): any {
    if (results.length === 0) {
      return { diversity_note: 'No results to analyze' };
    }

    const molecularWeights = results
      .map(r => r.properties.molecular_weight)
      .filter(mw => mw !== null && mw !== undefined);

    const alogpValues = results
      .map(r => r.properties.alogp)
      .filter(alogp => alogp !== null && alogp !== undefined);

    const moleculeTypes = results.reduce((acc: any, r: any) => {
      acc[r.molecule_type] = (acc[r.molecule_type] || 0) + 1;
      return acc;
    }, {});

    const maxPhases = results.reduce((acc: any, r: any) => {
      if (r.max_phase !== null && r.max_phase !== undefined) {
        acc[r.max_phase] = (acc[r.max_phase] || 0) + 1;
      }
      return acc;
    }, {});

    return {
      total_compounds: results.length,
      molecular_weight_range: molecularWeights.length > 0 ? {
        min: Math.min(...molecularWeights),
        max: Math.max(...molecularWeights),
        average: molecularWeights.reduce((a, b) => a + b, 0) / molecularWeights.length,
      } : null,
      alogp_range: alogpValues.length > 0 ? {
        min: Math.min(...alogpValues),
        max: Math.max(...alogpValues),
        average: alogpValues.reduce((a, b) => a + b, 0) / alogpValues.length,
      } : null,
      molecule_type_distribution: moleculeTypes,
      development_phase_distribution: maxPhases,
    };
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
    if (!args || typeof args.chembl_id !== 'string') {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid arguments: chembl_id is required');
    }

    try {
      // Get the main compound data first
      const compoundResponse = await this.apiClient.get(`/molecule/${args.chembl_id}.json`);
      const compound = compoundResponse.data;

      // Get cross-references from multiple ChEMBL endpoints
      const crossReferences = await this.gatherCrossReferences(args.chembl_id);

      // Organize references by database type
      const organizedReferences = this.organizeCrossReferences(crossReferences);

      const result = {
        chembl_id: compound.molecule_chembl_id,
        pref_name: compound.pref_name,
        molecule_type: compound.molecule_type,
        
        // Structural identifiers
        structural_identifiers: {
          canonical_smiles: compound.molecule_structures?.canonical_smiles,
          standard_inchi: compound.molecule_structures?.standard_inchi,
          standard_inchi_key: compound.molecule_structures?.standard_inchi_key,
        },

        // External database references organized by category
        external_references: organizedReferences,

        // Additional identifiers and aliases
        identifiers: {
          synonyms: compound.molecule_synonyms?.map((syn: any) => ({
            synonym: syn.molecule_synonym,
            syn_type: syn.syn_type,
          })) || [],
          trade_names: compound.molecule_synonyms?.filter((syn: any) => 
            syn.syn_type === 'TRADE_NAME'
          ).map((syn: any) => syn.molecule_synonym) || [],
        },

        // Hierarchy and classification
        hierarchy: compound.molecule_hierarchy || {},

        // Summary statistics
        reference_summary: {
          total_references: this.countTotalReferences(organizedReferences),
          database_count: Object.keys(organizedReferences).length,
          has_pdb_structures: organizedReferences.structural?.pdb?.length > 0,
          has_pubmed_refs: organizedReferences.literature?.pubmed?.length > 0,
          has_patent_refs: organizedReferences.patents?.length > 0,
        },

        // Usage recommendations
        usage_recommendations: this.generateUsageRecommendations(organizedReferences),
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get external references: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async gatherCrossReferences(chemblId: string): Promise<any[]> {
    const allReferences = [];

    try {
      // Get molecule cross-references
      const molXrefResponse = await this.apiClient.get('/molecule_xref.json', {
        params: {
          molecule_chembl_id: chemblId,
          limit: 100,
        },
      });
      allReferences.push(...(molXrefResponse.data.molecule_xrefs || []));
    } catch (error) {
      // Continue if this fails
    }

    try {
      // Get compound records which may have additional references
      const compoundResponse = await this.apiClient.get('/compound_record.json', {
        params: {
          molecule_chembl_id: chemblId,
          limit: 50,
        },
      });
      
      // Extract references from compound records
      const records = compoundResponse.data.compound_records || [];
      records.forEach((record: any) => {
        if (record.src_id) {
          allReferences.push({
            xref_src_db: record.src_description || 'Source Database',
            xref_id: record.src_compound_id,
            xref_name: record.compound_name,
          });
        }
      });
    } catch (error) {
      // Continue if this fails
    }

    return allReferences;
  }

  private organizeCrossReferences(references: any[]): any {
    const organized: any = {
      chemical_databases: {},
      structural: {},
      literature: {},
      patents: [],
      biological: {},
      other: [],
    };

    references.forEach((ref: any) => {
      const dbName = ref.xref_src_db?.toLowerCase() || '';
      const reference = {
        database: ref.xref_src_db,
        id: ref.xref_id,
        name: ref.xref_name,
        url: this.generateReferenceUrl(ref.xref_src_db, ref.xref_id),
      };

      // Categorize by database type
      if (dbName.includes('pubchem') || dbName.includes('chemspider') || dbName.includes('chebi')) {
        organized.chemical_databases[dbName] = organized.chemical_databases[dbName] || [];
        organized.chemical_databases[dbName].push(reference);
      } else if (dbName.includes('pdb') || dbName.includes('rcsb')) {
        organized.structural.pdb = organized.structural.pdb || [];
        organized.structural.pdb.push(reference);
      } else if (dbName.includes('pubmed') || dbName.includes('doi')) {
        organized.literature.pubmed = organized.literature.pubmed || [];
        organized.literature.pubmed.push(reference);
      } else if (dbName.includes('patent')) {
        organized.patents.push(reference);
      } else if (dbName.includes('uniprot') || dbName.includes('gene') || dbName.includes('kegg')) {
        organized.biological[dbName] = organized.biological[dbName] || [];
        organized.biological[dbName].push(reference);
      } else {
        organized.other.push(reference);
      }
    });

    return organized;
  }

  private generateReferenceUrl(database: string, id: string): string | null {
    if (!database || !id) return null;

    const dbLower = database.toLowerCase();
    
    // Common database URL patterns
    const urlMappings: { [key: string]: string } = {
      'pubchem': `https://pubchem.ncbi.nlm.nih.gov/compound/${id}`,
      'chemspider': `https://www.chemspider.com/Chemical-Structure.${id}.html`,
      'chebi': `https://www.ebi.ac.uk/chebi/searchId.do?chebiId=${id}`,
      'pdb': `https://www.rcsb.org/structure/${id}`,
      'pubmed': `https://pubmed.ncbi.nlm.nih.gov/${id}`,
      'uniprot': `https://www.uniprot.org/uniprot/${id}`,
      'kegg': `https://www.genome.jp/entry/${id}`,
    };

    for (const [key, urlPattern] of Object.entries(urlMappings)) {
      if (dbLower.includes(key)) {
        return urlPattern;
      }
    }

    return null;
  }

  private countTotalReferences(organized: any): number {
    let count = 0;
    
    Object.values(organized.chemical_databases).forEach((refs: any) => {
      count += Array.isArray(refs) ? refs.length : 0;
    });
    
    Object.values(organized.structural).forEach((refs: any) => {
      count += Array.isArray(refs) ? refs.length : 0;
    });
    
    Object.values(organized.literature).forEach((refs: any) => {
      count += Array.isArray(refs) ? refs.length : 0;
    });
    
    count += organized.patents?.length || 0;
    
    Object.values(organized.biological).forEach((refs: any) => {
      count += Array.isArray(refs) ? refs.length : 0;
    });
    
    count += organized.other?.length || 0;
    
    return count;
  }

  private generateUsageRecommendations(organized: any): string[] {
    const recommendations = [];

    if (organized.structural?.pdb?.length > 0) {
      recommendations.push('PDB structures available - suitable for structure-based drug design');
    }

    if (organized.chemical_databases?.pubchem?.length > 0) {
      recommendations.push('PubChem data available - check for additional bioactivity data');
    }

    if (organized.literature?.pubmed?.length > 0) {
      recommendations.push('Literature references available - review for mechanism and SAR data');
    }

    if (organized.biological?.uniprot?.length > 0) {
      recommendations.push('UniProt references available - check target protein information');
    }

    if (organized.patents?.length > 0) {
      recommendations.push('Patent information available - review for IP considerations');
    }

    if (recommendations.length === 0) {
      recommendations.push('Limited external references - consider expanding search to related compounds');
    }

    return recommendations;
  }

  private async handleAdvancedSearch(args: any) {
    if (!isValidPropertyFilterArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid advanced search arguments');
    }

    try {
      // Build query parameters for molecular property filtering
      const params: any = {
        limit: args.limit || 25,
      };

      // Add molecular weight filters
      if (args.min_mw !== undefined) {
        params['molecule_properties__mw_freebase__gte'] = args.min_mw;
      }
      if (args.max_mw !== undefined) {
        params['molecule_properties__mw_freebase__lte'] = args.max_mw;
      }

      // Add LogP filters
      if (args.min_logp !== undefined) {
        params['molecule_properties__alogp__gte'] = args.min_logp;
      }
      if (args.max_logp !== undefined) {
        params['molecule_properties__alogp__lte'] = args.max_logp;
      }

      // Add hydrogen bond donor/acceptor filters
      if (args.max_hbd !== undefined) {
        params['molecule_properties__hbd__lte'] = args.max_hbd;
      }
      if (args.max_hba !== undefined) {
        params['molecule_properties__hba__lte'] = args.max_hba;
      }

      // Execute the search
      const response = await this.apiClient.get('/molecule.json', { params });
      const compounds = response.data.molecules || [];

      // Analyze the results
      const analysis = this.analyzeAdvancedSearchResults(compounds, args);

      const result = {
        search_parameters: {
          molecular_weight_range: {
            min: args.min_mw || 'not specified',
            max: args.max_mw || 'not specified',
          },
          logp_range: {
            min: args.min_logp || 'not specified',
            max: args.max_logp || 'not specified',
          },
          hydrogen_bonding_limits: {
            max_donors: args.max_hbd || 'not specified',
            max_acceptors: args.max_hba || 'not specified',
          },
          limit: args.limit || 25,
        },
        results_summary: {
          total_found: compounds.length,
          compounds_analyzed: Math.min(compounds.length, args.limit || 25),
        },
        property_analysis: analysis,
        compounds: compounds.slice(0, args.limit || 25).map((compound: any) => ({
          chembl_id: compound.molecule_chembl_id,
          pref_name: compound.pref_name,
          molecule_type: compound.molecule_type,
          molecular_properties: {
            molecular_weight: compound.molecule_properties?.mw_freebase,
            alogp: compound.molecule_properties?.alogp,
            hbd: compound.molecule_properties?.hbd,
            hba: compound.molecule_properties?.hba,
            psa: compound.molecule_properties?.psa,
            rotatable_bonds: compound.molecule_properties?.rtb,
            lipinski_violations: compound.molecule_properties?.num_ro5_violations,
          },
          structures: {
            canonical_smiles: compound.molecule_structures?.canonical_smiles,
            standard_inchi_key: compound.molecule_structures?.standard_inchi_key,
          },
        })),
        search_insights: this.generateAdvancedSearchInsights(compounds, args),
        recommendations: this.generateAdvancedSearchRecommendations(compounds, args),
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(result, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to perform advanced search: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private analyzeAdvancedSearchResults(compounds: any[], searchParams: any): any {
    if (compounds.length === 0) {
      return {
        message: 'No compounds found matching the specified criteria',
        suggestions: 'Consider relaxing search parameters',
      };
    }

    const properties = compounds.map(c => c.molecule_properties).filter(p => p);
    
    const analysis = {
      molecular_weight: this.calculatePropertyStats(properties, 'mw_freebase'),
      alogp: this.calculatePropertyStats(properties, 'alogp'),
      hydrogen_bond_donors: this.calculatePropertyStats(properties, 'hbd'),
      hydrogen_bond_acceptors: this.calculatePropertyStats(properties, 'hba'),
      polar_surface_area: this.calculatePropertyStats(properties, 'psa'),
      rotatable_bonds: this.calculatePropertyStats(properties, 'rtb'),
      drug_likeness: {
        lipinski_compliant: properties.filter(p => (p.num_ro5_violations || 0) === 0).length,
        total_compounds: properties.length,
        compliance_rate: properties.length > 0 ? 
          (properties.filter(p => (p.num_ro5_violations || 0) === 0).length / properties.length * 100).toFixed(1) + '%' : '0%',
      },
    };

    return analysis;
  }

  private calculatePropertyStats(properties: any[], propertyName: string): any {
    const values = properties
      .map(p => p[propertyName])
      .filter(v => v !== null && v !== undefined && !isNaN(v))
      .map(v => parseFloat(v));

    if (values.length === 0) {
      return { message: 'No valid data available' };
    }

    values.sort((a, b) => a - b);
    
    return {
      count: values.length,
      min: values[0],
      max: values[values.length - 1],
      mean: (values.reduce((a, b) => a + b, 0) / values.length).toFixed(2),
      median: values.length % 2 === 0 
        ? ((values[values.length / 2 - 1] + values[values.length / 2]) / 2).toFixed(2)
        : values[Math.floor(values.length / 2)].toFixed(2),
    };
  }

  private generateAdvancedSearchInsights(compounds: any[], searchParams: any): string[] {
    const insights = [];

    if (compounds.length === 0) {
      insights.push('No compounds found - consider broadening search criteria');
      return insights;
    }

    const properties = compounds.map(c => c.molecule_properties).filter(p => p);
    
    // Molecular weight insights
    if (searchParams.min_mw || searchParams.max_mw) {
      const avgMw = properties.reduce((sum, p) => sum + (p.mw_freebase || 0), 0) / properties.length;
      insights.push(`Average molecular weight: ${avgMw.toFixed(1)} Da`);
    }

    // Drug-likeness insights
    const lipinskiCompliant = properties.filter(p => (p.num_ro5_violations || 0) === 0).length;
    const complianceRate = (lipinskiCompliant / properties.length * 100).toFixed(1);
    insights.push(`${complianceRate}% of compounds are Lipinski Rule of Five compliant`);

    // LogP distribution insights
    if (searchParams.min_logp || searchParams.max_logp) {
      const logpValues = properties.map(p => p.alogp).filter(v => v !== null && v !== undefined);
      if (logpValues.length > 0) {
        const avgLogp = logpValues.reduce((a, b) => a + b, 0) / logpValues.length;
        insights.push(`Average LogP: ${avgLogp.toFixed(2)} (lipophilicity indicator)`);
      }
    }

    // Diversity insights
    const moleculeTypes = [...new Set(compounds.map(c => c.molecule_type))];
    if (moleculeTypes.length > 1) {
      insights.push(`Chemical diversity: ${moleculeTypes.length} different molecule types found`);
    }

    return insights;
  }

  private generateAdvancedSearchRecommendations(compounds: any[], searchParams: any): string[] {
    const recommendations = [];

    if (compounds.length === 0) {
      recommendations.push('Try increasing molecular weight range');
      recommendations.push('Consider relaxing LogP constraints');
      recommendations.push('Remove hydrogen bonding limitations');
      return recommendations;
    }

    const properties = compounds.map(c => c.molecule_properties).filter(p => p);
    
    // Check if results are heavily skewed
    if (compounds.length < (searchParams.limit || 25) / 2) {
      recommendations.push('Consider broadening search parameters for more diverse results');
    }

    // Drug-likeness recommendations
    const lipinskiCompliant = properties.filter(p => (p.num_ro5_violations || 0) === 0).length;
    const complianceRate = lipinskiCompliant / properties.length;
    
    if (complianceRate < 0.5) {
      recommendations.push('Many compounds violate Lipinski rules - consider tightening MW/LogP constraints for drug-like properties');
    } else if (complianceRate > 0.9) {
      recommendations.push('High drug-likeness compliance - good for lead optimization studies');
    }

    // Specific property recommendations
    if (searchParams.min_mw && searchParams.max_mw) {
      const range = searchParams.max_mw - searchParams.min_mw;
      if (range < 100) {
        recommendations.push('Narrow MW range may limit chemical diversity - consider expanding');
      }
    }

    if (searchParams.max_hbd !== undefined && searchParams.max_hbd < 3) {
      recommendations.push('Low HBD limit may exclude important pharmacophores');
    }

    // Analysis recommendations
    recommendations.push('Use substructure_search for scaffold-based filtering');
    recommendations.push('Apply calculate_descriptors for detailed property analysis');
    recommendations.push('Consider get_external_references for additional compound information');

    return recommendations;
  }

  async run() {
    const transport = new StdioServerTransport();
    await this.server.connect(transport);
    console.error('ChEMBL MCP server running on stdio');
  }
}

const server = new ChEMBLServer();
server.run().catch(console.error);
