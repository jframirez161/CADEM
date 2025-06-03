# api/main.py
from fastapi import FastAPI, HTTPException, Body
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field # Field for validation
from typing import List, Dict, Any, Optional # Optional for potentially nullable fields


# Relative import assuming 'simulation' is a package in the same directory level as main.py
# and 'api' is the root for this part of the project.
from simulation.run_model import execute_single_simulation

app = FastAPI(title="CADEM Rumen Model API")

# --- CORS Middleware ---
# Adjust origins as needed for your frontend's URL
origins = [
    "https://cadem-frontend.onrender.com", # Your deployed frontend
    "http://localhost:5173",              # Your local Vite dev server (common default)
    "http://localhost:3000",              # Another common local dev server port
    # Add any other origins you need to allow (e.g., staging environment)
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["POST", "GET"], # Specify allowed methods
    allow_headers=["*"],
)

# --- Pydantic Models for Request/Response Validation ---
class SimulationParamsAPI(BaseModel): # Renamed to avoid conflict with other uses
    hours: int = Field(default=24, ge=1, description="Total simulation duration in hours")
    points_per_hour: int = Field(default=1000, ge=1, le=5000, description="Number of data points to output per hour")

class DietParamsAPI(BaseModel): # Renamed
    DMI: float = Field(..., gt=0, description="Dry Matter Intake in kg/day")
    NDF: float = Field(..., ge=0, description="Neutral Detergent Fiber in g/kg DM")
    St: float = Field(..., ge=0, description="Starch in g/kg DM")
    WSC: float = Field(..., ge=0, description="Water Soluble Carbohydrates in g/kg DM")
    Acin: float = Field(..., ge=0, description="Acetate input in g/kg DM")
    Prin: float = Field(..., ge=0, description="Propionate input in g/kg DM")
    Buin: float = Field(..., ge=0, description="Butyrate input in g/kg DM")
    Lain: float = Field(..., ge=0, description="Lactate input in g/kg DM")
    # Add other essential diet nutrients the model expects, e.g.:
    CP: Optional[float] = Field(None, ge=0, description="Crude Protein in g/kg DM (optional)")
    Fat: Optional[float] = Field(None, ge=0, description="Fat (Ether Extract) in g/kg DM (optional)")
    Ash: Optional[float] = Field(None, ge=0, description="Ash in g/kg DM (optional)")
    # Make them optional if the frontend might not always send them,
    # and handle defaults or errors in get_model_constants if they are truly required.

class SimulationRequest(BaseModel):
    simulation_params: SimulationParamsAPI
    diet_params: DietParamsAPI
    # Example for optional complex input:
    # feed_intake_pattern: Optional[List[Dict[str, float]]] = None

@app.post("/api/simulate", summary="Run Rumen Simulation")
async def simulate_rumen_model_endpoint(request_data: SimulationRequest = Body(...)):
    """
    Accepts simulation and diet parameters, runs the rumen model,
    and returns the time-series results.
    """
    print(f"Received API simulation request: {request_data.dict(exclude_none=True)}")

    # Pydantic models automatically convert to dicts when calling .dict()
    result = execute_single_simulation(
        sim_params=request_data.simulation_params.dict(),
        diet_params=request_data.diet_params.dict(exclude_none=True) # Exclude Nones for cleaner diet_params
    )

    if not result.get("success"):
        # Log the error on the server for more details if needed
        print(f"API call to /api/simulate failed: {result.get('error', 'Unknown simulation error')}")
        raise HTTPException(
            status_code=500, 
            detail=result.get("error", "Simulation failed internally.")
        )
    
    return {
        "success": True,
        "results": result.get("results", []),
        "message": "Simulation successful."
    }

@app.get("/api/health", summary="API Health Check")
async def health_check():
    """
    Simple health check endpoint.
    """
    return {"status": "CADEM Rumen Model API is running"}

if __name__ == "__main__":
    import uvicorn
    # For local development: uvicorn main:app --reload --host 0.0.0.0 --port 8000
    # The host "0.0.0.0" makes it accessible from other devices on the network.
    # The port 8000 is a common choice for APIs.
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)